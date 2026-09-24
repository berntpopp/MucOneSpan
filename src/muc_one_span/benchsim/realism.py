"""Realism metrics of simulated reads, measured against their per-read truth.

`read_metrics` measures one case; `aggregate` pools several matched cases;
`realism_targets.compare` scores either against the public PRJEB92208 targets
(spec section 4 tolerances). Measurements mirror the real-data profiler that
produced the targets (``realprofile``) so both sides use the same definitions:

- Orientation: ``-`` reads are reverse-complemented into haplotype orientation
  before alignment, so ``+``/``-`` match the targets' ``_meta.strand_convention``
  (``+`` = C-runs read as C). Each ``full``/``fragment`` read is aligned
  globally (edlib ``mode="NW"``, ``task="path"``) to its truth source interval.
  Source frames: genomic truth coordinates index the haplotype; amplicon
  coordinates index the primer-to-primer amplicon (`amplicon_sources`).
- Error rates per strand (``+``, ``-``, ``all``): mismatch, inserted and
  deleted bases per source (reference) base, and their sum ``error_rate``.
  Targets: ``error_rates_per_ref_base["<strand>:<mismatch|ins|del|total>"]``.
- C7 correct-length fraction per strand: every run of exactly 7 C in the
  source (not touching an interval end) is measured through the alignment
  (read bases between the run's aligned bounds, extended by up to 3 adjacent
  read C) and counted correct when 7 are observed. Target:
  ``hp_P_obs_given_true["C7|<strand>"]["p_correct"]``.
- Amplicon products (``full``/``smear``/``chimera``/``concatemer``): the
  offset ``d`` = read length minus the nearest allele's amplicon length
  (per-haplotype ``full`` source length). ``span_off_gt1unit_frac`` (|d| > 1.5
  units), ``span_between_alleles_frac`` and ``span_below_short_frac`` are the
  observable smear products, as in the targets' same-named keys; the 15 bp
  offset histogram (bins ``(d + 7) // 15`` clamped to [-12, 4], split by
  allele size < / >= 55 units) matches ``span_offset_pmf_15bp_bins``.
  ``smear_frac``/``chimera_frac`` (truth kinds among products) are reported
  only; ``offtarget_frac`` (of all reads) is scored against
  ``category_frac["off_target"]``.
- Allele ratio: on-peak products (|d| <= max(30 bp, 1.2 %)) per allele give
  ``log_ratio`` = ln(n_long / n_short) over ``delta_units``; the per-case
  ``log_ratio_slope`` is their quotient and `aggregate` fits a through-origin
  slope across cases (target ``allelic_ratio["through_origin_b_per_unit"]``).
- ``spanning_frac`` (genomic, when ``vntr_span`` is given): fragments covering
  the haplotype's VNTR interval over all reads (``category_frac["spanning"]``).

Known differences from the real profiler: it measured errors only between the
motif-1/motif-9 anchors of on-peak spanning reads, whereas the simulation uses
the whole source interval of every ``full``/``fragment`` read; amplicon lengths
here include primer-to-anchor flanks, which shifts the 55-unit size split by
about two units.
"""

from __future__ import annotations

import importlib
import math
import re
import statistics
from collections import Counter
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

from muc_one_span.evaluation.truth import fasta_records
from muc_one_span.nomenclature import revcomp

from .read_truth import ReadTruth

UNIT = 60
STRANDS = ("+", "-")
ALIGNED_KINDS = ("full", "fragment")
PRODUCT_KINDS = ("full", "smear", "chimera", "concatemer")
RATE_KEYS = {
    "mismatch_rate": ("mismatch",),
    "insertion_rate": ("ins",),
    "deletion_rate": ("del",),
    "error_rate": ("mismatch", "ins", "del"),
}
FRACTION_KEYS = (
    "smear_frac",
    "chimera_frac",
    "offtarget_frac",
    "span_off_gt1unit_frac",
    "span_between_alleles_frac",
    "span_below_short_frac",
    "spanning_frac",
)
BIN_MIN, BIN_MAX = -12, 4
SIZE_SPLIT_UNITS = 55
C7 = 7
COUNT_KEYS = ("ref_bases", "mismatch", "ins", "del")
_CIGAR = re.compile(r"(\d+)([=XID])")


def _edlib() -> ModuleType:
    try:
        return importlib.import_module("edlib")
    except ImportError as exc:
        raise ImportError(
            "realism metrics need edlib; install the optional extra "
            "'muc-one-span[bench]' (uv sync --extra bench)"
        ) from exc


def _fastq(path: Path) -> Iterator[tuple[str, str]]:
    with path.open() as handle:
        while header := handle.readline():
            seq = handle.readline().strip().upper()
            handle.readline()
            handle.readline()
            yield header[1:].split(maxsplit=1)[0], seq


def _align(edlib: ModuleType, read: str, ref: str) -> tuple[Counter[str], list[int]]:
    """Global alignment error counts and a source-index -> read-index map."""
    counts: Counter[str] = Counter(ref_bases=len(ref))
    ref2read = [0] * (len(ref) + 1)
    qi = ri = 0
    cigar = edlib.align(read, ref, mode="NW", task="path")["cigar"] or ""
    for num, op in _CIGAR.findall(cigar):
        n = int(num)
        if op == "I":
            counts["ins"] += n
            qi += n
            continue
        for k in range(n):
            ref2read[ri + k] = qi + (k if op != "D" else 0)
        ri += n
        if op == "D":
            counts["del"] += n
        else:
            qi += n
            counts["mismatch"] += n if op == "X" else 0
    ref2read[len(ref)] = qi
    return counts, ref2read


def _c7_calls(read: str, ref: str, ref2read: list[int]) -> list[bool]:
    """Per C7 source run: whether the read shows exactly 7 C there."""
    calls = []
    for match in re.finditer(r"C+", ref):
        s, e = match.span()
        if e - s != C7 or s == 0 or e >= len(ref):
            continue
        rs, re_ = ref2read[s], ref2read[e]
        obs = read[rs:re_].count("C")
        if obs:
            i = rs - 1
            while i >= max(rs - 3, 0) and read[i] == "C":
                obs, i = obs + 1, i - 1
            i = re_
            while i < min(re_ + 3, len(read)) and read[i] == "C":
                obs, i = obs + 1, i + 1
        calls.append(obs == C7)
    return calls


def _count_table(counts: dict[str, Counter[str]]) -> dict[str, dict[str, int]]:
    return {s: {k: c[k] for k in COUNT_KEYS} for s, c in counts.items()}


def _ratio(num: float, den: float) -> float | None:
    return num / den if den else None


def _derived(error_counts: dict[str, Counter[str]], c7: dict[str, list[int]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, parts in RATE_KEYS.items():
        out[key] = {
            s: _ratio(sum(error_counts[s][p] for p in parts), error_counts[s]["ref_bases"])
            for s in (*STRANDS, "all")
        }
    out["c7_correct"] = {s: _ratio(ok, n) for s, (ok, n) in c7.items()}
    return out


def _alleles(rows: Sequence[ReadTruth]) -> list[int]:
    per_hap: dict[int, Counter[int]] = {}
    for r in rows:
        if r.kind == "full":
            per_hap.setdefault(r.hap, Counter())[r.src_end - r.src_start] += 1
    return sorted(c.most_common(1)[0][0] for c in per_hap.values())


def _amplicon(rows: Sequence[ReadTruth], lengths: dict[str, int]) -> dict[str, Any]:
    kinds = Counter(r.kind for r in rows)
    products = [r for r in rows if r.kind in PRODUCT_KINDS]
    if not (products or kinds["offtarget"]):
        return {}  # genomic fragments: no amplicon products
    alleles = _alleles(rows)
    if not alleles:
        return {"offtarget_frac": _ratio(kinds["offtarget"], len(rows))}
    lo, hi = alleles[0], alleles[-1]
    hist = {"lt55u": [0] * (BIN_MAX - BIN_MIN + 1), "ge55u": [0] * (BIN_MAX - BIN_MIN + 1)}
    on_peak: Counter[int] = Counter()
    off = between = below = 0
    for r in products:
        span = lengths[r.read_id]
        allele = min(alleles, key=lambda a: abs(span - a))
        d = span - allele
        key = "lt55u" if allele / UNIT < SIZE_SPLIT_UNITS else "ge55u"
        hist[key][max(BIN_MIN, min(BIN_MAX, (d + 7) // 15)) - BIN_MIN] += 1
        on_peak[allele] += abs(d) <= max(30.0, 0.012 * allele)
        off += abs(d) > 1.5 * UNIT
        between += lo + 1.5 * UNIT < span < hi - 1.5 * UNIT
        below += span < lo - 1.5 * UNIT
    n = len(products)
    out: dict[str, Any] = {
        "smear_frac": _ratio(kinds["smear"], n),
        "chimera_frac": _ratio(kinds["chimera"], n),
        "offtarget_frac": _ratio(kinds["offtarget"], len(rows)),
        "span_off_gt1unit_frac": _ratio(off, n),
        "span_between_alleles_frac": _ratio(between, n),
        "span_below_short_frac": _ratio(below, n),
        "span_offset_hist": hist,
    }
    delta = (hi - lo) / UNIT
    if delta >= 1 and on_peak[lo] and on_peak[hi]:
        log_ratio = math.log(on_peak[hi] / on_peak[lo])
        out["allele_ratio"] = {"delta_units": delta, "log_ratio": log_ratio}
        out["log_ratio_slope"] = log_ratio / delta
    return out


def _empty() -> dict[str, Any]:
    return dict.fromkeys(FRACTION_KEYS) | {
        "allele_ratio": None,
        "log_ratio_slope": None,
        "span_offset_hist": None,
    }


def haplotype_sources(truth_fasta: Path) -> dict[int, str]:
    """Haplotype id -> sequence from a MucOneUp truth FASTA (``haplotype_<n>``)."""
    records = fasta_records(truth_fasta)
    out = {}
    for index, (name, seq) in enumerate(records.items(), 1):
        match = re.search(r"haplotype_(\d+)", name)
        out[int(match.group(1)) if match else index] = seq.upper()
    return out


def amplicon_sources(
    haplotypes: dict[int, str], forward_primer: str, reverse_primer: str
) -> dict[int, str]:
    """Primer-to-primer amplicon per haplotype (MucOneUp amplicon source frame).

    Amplicon read truth coordinates are relative to the extracted amplicon
    (forward primer start to reverse primer end, both included), which is
    reproduced here with the same exact, unique primer matching.

    Raises:
        ValueError: If a primer site is missing or not unique on a haplotype.
    """
    fwd, rev = forward_primer.upper(), revcomp(reverse_primer.upper())
    out = {}
    for hap, seq in haplotypes.items():
        f_sites = [m.start() for m in re.finditer(f"(?={fwd})", seq)]
        r_sites = [m.start() for m in re.finditer(f"(?={rev})", seq)]
        if len(f_sites) != 1 or len(r_sites) != 1 or r_sites[0] < f_sites[0]:
            raise ValueError(f"haplotype {hap}: primer sites not unique ({f_sites}, {r_sites})")
        out[hap] = seq[f_sites[0] : r_sites[0] + len(rev)]
    return out


def read_metrics(
    fastq: Path,
    truth_rows: Sequence[ReadTruth],
    sources: dict[int, str],
    vntr_span: dict[int, tuple[int, int]] | None = None,
) -> dict[str, Any]:
    """Measure the realism metrics of one simulated case.

    Args:
        fastq: The case's reads FASTQ.
        truth_rows: Its per-read truth rows (`load_read_truth`).
        sources: Haplotype id -> sequence in the read truth frame: the
            haplotype for genomic fragments (`haplotype_sources`), its
            primer-to-primer amplicon for amplicon reads (`amplicon_sources`).
        vntr_span: Optional haplotype -> `[start, end)` VNTR interval; when
            given, `spanning_frac` counts fragments covering it (genomic).

    Returns:
        Per-case metrics: raw `error_counts`/`c7_counts` (for pooling), the
        rates and `c7_correct` per strand (`None` without data), amplicon
        fractions, `span_offset_hist`, `allele_ratio`, `log_ratio_slope` and
        `spanning_frac` (`None` where not applicable).

    Raises:
        ImportError: If edlib (extra `bench`) is not installed.
        ValueError: If a truth read is missing from the FASTQ or its
            haplotype has no source sequence.
    """
    edlib = _edlib()
    seqs = dict(_fastq(fastq))
    missing = [r.read_id for r in truth_rows if r.read_id not in seqs]
    if missing:
        raise ValueError(f"{len(missing)} truth reads are not in the FASTQ, e.g. {missing[0]}")
    counts: dict[str, Counter[str]] = {s: Counter() for s in (*STRANDS, "all")}
    c7 = {s: [0, 0] for s in (*STRANDS, "both")}
    n_aligned = dict.fromkeys(STRANDS, 0)
    for r in truth_rows:
        if r.kind not in ALIGNED_KINDS:
            continue
        if r.hap not in sources:
            raise ValueError(f"no source sequence for haplotype {r.hap} ({r.read_id})")
        read = seqs[r.read_id] if r.strand == "+" else revcomp(seqs[r.read_id])
        ref = sources[r.hap][r.src_start : r.src_end].upper()
        if not read or not ref:
            continue
        found, ref2read = _align(edlib, read, ref)
        counts[r.strand].update(found)
        counts["all"].update(found)
        n_aligned[r.strand] += 1
        for ok in _c7_calls(read, ref, ref2read):
            for key in (r.strand, "both"):
                c7[key][0] += ok
                c7[key][1] += 1
    out = _empty() | _amplicon(truth_rows, {k: len(v) for k, v in seqs.items()})
    if vntr_span is not None:
        spanning = sum(
            r.kind == "fragment"
            and r.hap in vntr_span
            and r.src_start <= vntr_span[r.hap][0]
            and r.src_end >= vntr_span[r.hap][1]
            for r in truth_rows
        )
        out["spanning_frac"] = _ratio(spanning, len(truth_rows))
    out |= _derived(counts, c7)
    out["error_counts"] = _count_table(counts)
    out["c7_counts"] = c7
    out["n_aligned"] = n_aligned
    out["n_reads"] = len(truth_rows)
    return out


def _spread(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    return {
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
        "n": len(values),
    }


def aggregate(cases: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Pool per-case `read_metrics` of matched cases into across-case metrics.

    Error and C7 counts and the offset histograms are summed; fractions become
    `{min, median, max, n}` over cases (scored as range + median by
    `compare`); `log_ratio_slope` is the through-origin least-squares slope of
    `log_ratio` on `delta_units` across cases with two resolvable alleles.

    Args:
        cases: Per-case metric dicts from `read_metrics`.

    Returns:
        An aggregate metrics dict accepted by `realism_targets.compare`.
    """
    counts: dict[str, Counter[str]] = {s: Counter() for s in (*STRANDS, "all")}
    c7 = {s: [0, 0] for s in (*STRANDS, "both")}
    hist = {k: [0] * (BIN_MAX - BIN_MIN + 1) for k in ("lt55u", "ge55u")}
    points = []
    for case in cases:
        for s, c in case["error_counts"].items():
            counts[s].update(c)
        for s, (ok, n) in case["c7_counts"].items():
            c7[s][0] += ok
            c7[s][1] += n
        for k, h in (case.get("span_offset_hist") or {}).items():
            hist[k] = [a + b for a, b in zip(hist[k], h, strict=True)]
        if case.get("allele_ratio"):
            points.append((case["allele_ratio"]["delta_units"], case["allele_ratio"]["log_ratio"]))
    out: dict[str, Any] = {"n_cases": len(cases)}
    for key in FRACTION_KEYS:
        out[key] = _spread([c[key] for c in cases if c.get(key) is not None])
    sxx = sum(x * x for x, _ in points)
    out["log_ratio_slope"] = _ratio(sum(x * y for x, y in points), sxx)
    out["n_allele_ratio_cases"] = len(points)
    out["span_offset_hist"] = hist
    out |= _derived(counts, c7)
    out["error_counts"] = _count_table(counts)
    out["c7_counts"] = c7
    return out
