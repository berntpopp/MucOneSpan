"""Realism metrics of simulated reads, measured on the VNTR against per-read truth.

`read_metrics` measures one case (`case_metrics` loads it from a case
directory); `aggregate` combines matched cases; `realism_targets.compare`
scores either against the public PRJEB92208 targets (spec section 4). The
definitions mirror the real-data profiler that produced the targets
(``realprofile``), which measured anchor-to-anchor (motif-1 start to motif-9
end) spans of spanning reads:

- Scope: only ``full``/``fragment`` reads whose source interval covers the
  haplotype's VNTR (`geometry`, read-truth source frame) are scored. Each is
  aligned globally (edlib ``mode="NW"``, ``task="path"``) to its source
  interval; the read segment aligned to the VNTR is then realigned to the VNTR
  alone, and errors and C7 are counted on that alignment.
- Orientation: ``-`` reads are reverse-complemented into haplotype orientation,
  so ``+``/``-`` match the targets' ``_meta.strand_convention`` (``+`` =
  C-runs read as C).
- Error rates per strand (``+``, ``-``, ``all``): mismatch, inserted and
  deleted bases per VNTR base, and their sum ``error_rate``. Per case these
  are pooled over reads; `aggregate` reports the median of per-case rates,
  as the target median is over per-library rates.
- C7 correct-length fraction per strand: every run of exactly
  ``c7_run_length`` [7] C inside the VNTR (not touching its ends) is read
  through the alignment (read bases between the run's aligned bounds, extended
  by up to ``c7_extend_max`` [3] adjacent read C) and is correct when the run
  length is observed. As in the profiler, a run with no C in its
  aligned window (obs 0) is not extended and counts as incorrect. Counts are
  pooled, as in the target's ``hp_P_obs_given_true``.
- Amplicon products (``full``/``smear``/``chimera``/``concatemer``) are the
  simulated spanning reads, the targets' denominator for ``span_*`` fractions.
  Their VNTR span is the aligned segment length for ``full`` reads and read
  length minus the haplotype's amplicon flanks otherwise (smear, chimera and
  concatemer junctions lie inside the VNTR). ``d`` = span minus the nearest
  allele's VNTR length gives ``span_off_gt1unit_frac`` (|d| >
  ``off_peak_units`` [1.5] units), ``span_between_alleles_frac``,
  ``span_below_short_frac`` and the offset histogram (bins ``(d + bp // 2) //
  bp`` for ``offset_bin_bp`` [15], clamped to [``offset_bin_min``,
  ``offset_bin_max``] [-12, 4], allele size < / >= ``size_split_units`` [55]
  units) for ``span_offset_pmf_15bp_bins``. ``smear_frac``/``chimera_frac``
  (truth kinds over products) are reported only; ``offtarget_frac`` is over
  all reads, like ``category_frac["off_target"]``.
- Allele ratio: on-peak products (|d| <= max(``on_peak_min_bp`` [30 bp],
  ``on_peak_rel`` [1.2 %] of the allele)) per allele give
  ``log_ratio`` = ln(n_long / n_short) over ``delta_units``; the per-case
  ``log_ratio_slope`` is their quotient and `aggregate` fits a through-origin
  slope across cases (target ``allelic_ratio["through_origin_b_per_unit"]``).
- ``spanning_frac`` (genomic): VNTR-spanning fragments over fragments that
  overlap the VNTR interval.

Remaining differences from the profiler: it used a subsample of on-peak reads
for errors and assigned reads to alleles by length, not by truth; its
``spanning`` category also counts reads outside the VNTR in its denominator
(the pulled region), so ``spanning_frac`` here is an upper-side comparison.
"""

from __future__ import annotations

import importlib
import json
import math
import re
import statistics
from collections import Counter
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation.truth import load_truth
from muc_one_span.nomenclature import revcomp

from .bench_config import DEFAULT_BENCH_CONFIG, RealismConfig
from .generate import FASTQ
from .geometry import (
    GENOMIC,
    case_geometry,
    flank_sequences,
    haplotype_sequences,
    primer_pair,
    source_inputs,
    vntr_bounds,
)
from .read_truth import ReadTruth, load_read_truth

UNIT = load_repeat_dictionary().repeat_length_bp  # repeat-unit length (bundled dictionary)
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
    """Global alignment error counts and a reference-index -> read-index map."""
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


def _c7_calls(read: str, ref: str, ref2read: list[int], cfg: RealismConfig) -> list[bool]:
    """Per C7 reference run: whether the read shows exactly ``c7_run_length`` C there."""
    run, extend = cfg.c7_run_length, cfg.c7_extend_max
    calls = []
    for match in re.finditer(r"C+", ref):
        s, e = match.span()
        if e - s != run or s == 0 or e >= len(ref):
            continue
        rs, re_ = ref2read[s], ref2read[e]
        obs = read[rs:re_].count("C")
        if obs:  # profiler parity: an empty window is not extended
            i = rs - 1
            while i >= max(rs - extend, 0) and read[i] == "C":
                obs, i = obs + 1, i - 1
            i = re_
            while i < min(re_ + extend, len(read)) and read[i] == "C":
                obs, i = obs + 1, i + 1
        calls.append(obs == run)
    return calls


def _count_table(counts: dict[str, Counter[str]]) -> dict[str, dict[str, int]]:
    return {s: {k: c[k] for k in COUNT_KEYS} for s, c in counts.items()}


def _ratio(num: float, den: float) -> float | None:
    return num / den if den else None


def _rates(error_counts: dict[str, Counter[str]]) -> dict[str, Any]:
    return {
        key: {
            s: _ratio(sum(error_counts[s][p] for p in parts), error_counts[s]["ref_bases"])
            for s in (*STRANDS, "all")
        }
        for key, parts in RATE_KEYS.items()
    }


def _amplicon(
    rows: Sequence[ReadTruth],
    lengths: dict[str, int],
    exact: dict[str, int],
    sources: dict[int, str],
    vntr: dict[int, tuple[int, int]],
    cfg: RealismConfig,
) -> dict[str, Any]:
    kinds = Counter(r.kind for r in rows)
    products = [r for r in rows if r.kind in PRODUCT_KINDS]
    if not (products or kinds["offtarget"]):
        return {}  # genomic fragments: no amplicon products
    allele = {h: e - s for h, (s, e) in vntr.items()}
    flanks = {h: len(sources[h]) - n for h, n in allele.items() if h in sources}
    alleles = sorted(set(allele.values()))
    lo, hi = alleles[0], alleles[-1]
    small, large = cfg.size_keys
    hist = {small: [0] * cfg.n_bins, large: [0] * cfg.n_bins}
    bp, off_bp = cfg.offset_bin_bp, cfg.off_peak_units * UNIT
    on_peak: Counter[int] = Counter()
    off = between = below = 0
    for r in products:
        span = exact.get(r.read_id, lengths[r.read_id] - flanks.get(r.hap, 0))
        near = min(alleles, key=lambda a: abs(span - a))
        d = span - near
        key = small if near / UNIT < cfg.size_split_units else large
        b = max(cfg.offset_bin_min, min(cfg.offset_bin_max, (d + bp // 2) // bp))
        hist[key][b - cfg.offset_bin_min] += 1
        on_peak[near] += abs(d) <= max(cfg.on_peak_min_bp, cfg.on_peak_rel * near)
        off += abs(d) > off_bp
        between += lo + off_bp < span < hi - off_bp
        below += span < lo - off_bp
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


def read_metrics(
    fastq: Path,
    truth_rows: Sequence[ReadTruth],
    sources: dict[int, str],
    vntr: dict[int, tuple[int, int]],
    config: RealismConfig = DEFAULT_BENCH_CONFIG.realism,
) -> dict[str, Any]:
    """Measure the realism metrics of one simulated case on the VNTR interval.

    Args:
        fastq: The case's reads FASTQ.
        truth_rows: Its per-read truth rows (`load_read_truth`).
        sources: Haplotype id -> sequence in the read-truth source frame
            (`geometry.source_inputs`).
        vntr: Haplotype id -> VNTR `[start, end)` in the same frame.
        config: Metric definitions (C7 run, histogram bins, peak widths).

    Returns:
        Per-case metrics: raw `error_counts`/`c7_counts`, the rates and
        `c7_correct` per strand (`None` without data), amplicon fractions,
        `span_offset_hist`, `allele_ratio`, `log_ratio_slope`,
        `spanning_frac`, `n_aligned` (scored spanning reads) and `n_reads`.

    Raises:
        ImportError: If edlib (extra `bench`) is not installed.
        ValueError: If a truth read is missing from the FASTQ or its
            haplotype has no source sequence or VNTR bounds.
    """
    edlib = _edlib()
    seqs = dict(_fastq(fastq))
    missing = [r.read_id for r in truth_rows if r.read_id not in seqs]
    if missing:
        raise ValueError(f"{len(missing)} truth reads are not in the FASTQ, e.g. {missing[0]}")
    counts: dict[str, Counter[str]] = {s: Counter() for s in (*STRANDS, "all")}
    c7 = {s: [0, 0] for s in (*STRANDS, "both")}
    n_aligned = dict.fromkeys(STRANDS, 0)
    exact: dict[str, int] = {}
    overlap = spanning = 0
    for r in truth_rows:
        if r.kind not in ALIGNED_KINDS:
            continue
        if r.hap not in sources or r.hap not in vntr:
            raise ValueError(f"no source sequence or VNTR for haplotype {r.hap} ({r.read_id})")
        lo, hi = vntr[r.hap]
        overlap += r.src_start < hi and r.src_end > lo
        if not (r.src_start <= lo and r.src_end >= hi):
            continue
        spanning += 1
        read = seqs[r.read_id] if r.strand == "+" else revcomp(seqs[r.read_id])
        ref = sources[r.hap][r.src_start : r.src_end].upper()
        _, ref2read = _align(edlib, read, ref)
        segment = read[ref2read[lo - r.src_start] : ref2read[hi - r.src_start]]
        exact[r.read_id] = len(segment)
        vntr_ref = ref[lo - r.src_start : hi - r.src_start]
        found, seg2read = _align(edlib, segment, vntr_ref)
        counts[r.strand].update(found)
        counts["all"].update(found)
        n_aligned[r.strand] += 1
        for ok in _c7_calls(segment, vntr_ref, seg2read, config):
            for key in (r.strand, "both"):
                c7[key][0] += ok
                c7[key][1] += 1
    lengths = {k: len(v) for k, v in seqs.items()}
    out = _empty() | _amplicon(truth_rows, lengths, exact, sources, vntr, config)
    if any(r.kind == "fragment" for r in truth_rows):
        out["spanning_frac"] = _ratio(spanning, overlap)
    out |= _rates(counts)
    out["c7_correct"] = {s: _ratio(ok, n) for s, (ok, n) in c7.items()}
    out["error_counts"] = _count_table(counts)
    out["c7_counts"] = c7
    out["n_aligned"] = n_aligned
    out["n_reads"] = len(truth_rows)
    return out


def _single(root: Path, pattern: str) -> Path:
    found = sorted(root.glob(pattern))
    if len(found) != 1:
        raise ValueError(f"expected one {pattern} in {root}, found {len(found)}")
    return found[0]


def case_metrics(
    case_dir: Path,
    muconeup_config: Path | None = None,
    flank_fasta: Path | None = None,
    config: RealismConfig = DEFAULT_BENCH_CONFIG.realism,
) -> dict[str, Any]:
    """`read_metrics` for one generated case directory.

    Geometry comes from ``case.json["geometry"]``; cases generated before it
    was recorded are rebuilt from the truth and `muconeup_config` (primers).

    Args:
        case_dir: ``<out_root>/<split>/<design_id>``.
        muconeup_config: MucOneUp config, needed only for older amplicon cases.
        flank_fasta: The ``--flank-fasta`` used at generation (genomic).
        config: Metric definitions passed to `read_metrics`.

    Raises:
        ValueError: If an older amplicon case has no `muconeup_config`, or
            the flank lengths differ from the recorded ones.
    """
    case = json.loads((case_dir / "case.json").read_text())
    truth_dir = case_dir / "truth"
    haplotypes = haplotype_sequences(_single(truth_dir, "*.simulated.fa"))
    flanks = flank_sequences(flank_fasta)
    geometry = case.get("geometry")
    if geometry is None:
        genomic = case["profile"] == GENOMIC
        if not genomic and muconeup_config is None:
            raise ValueError(f"{case_dir}: case.json has no geometry; pass muconeup_config")
        rd = load_repeat_dictionary()
        geometry = case_geometry(
            case["profile"],
            haplotypes,
            vntr_bounds(load_truth(truth_dir, rd), rd.flanking_left, haplotypes),
            primers=None if muconeup_config is None else primer_pair(muconeup_config),
            flank_ext=(len(flanks[0]), len(flanks[1])),
        )
    sources, vntr = source_inputs(geometry, haplotypes, flanks)
    rows = load_read_truth(_single(truth_dir, "*_read_truth.tsv.gz"))
    fastq = case_dir / "reads" / FASTQ[case["profile"]].format(case["design_id"])
    return read_metrics(fastq, rows, sources, vntr, config)


def _spread(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    return {
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
        "n": len(values),
    }


def _median(values: list[float | None]) -> float | None:
    present = [v for v in values if v is not None]
    return statistics.median(present) if present else None


def aggregate(
    cases: Sequence[dict[str, Any]], config: RealismConfig = DEFAULT_BENCH_CONFIG.realism
) -> dict[str, Any]:
    """Combine per-case `read_metrics` of matched cases into across-case metrics.

    Error rates are the median of per-case rates (the target median is over
    per-library rates); C7 and error counts and the offset histograms are
    pooled; fractions become `{min, median, max, n}` over cases (range +
    median scoring in `compare`); `log_ratio_slope` is the through-origin
    least-squares slope of `log_ratio` on `delta_units` over cases.

    Args:
        cases: Per-case metric dicts from `read_metrics`.
        config: The histogram layout the cases were measured with.

    Returns:
        An aggregate metrics dict accepted by `realism_targets.compare`.
    """
    counts: dict[str, Counter[str]] = {s: Counter() for s in (*STRANDS, "all")}
    c7 = {s: [0, 0] for s in (*STRANDS, "both")}
    hist = {k: [0] * config.n_bins for k in config.size_keys}
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
    for key in RATE_KEYS:
        out[key] = {s: _median([c[key][s] for c in cases]) for s in (*STRANDS, "all")}
    for key in FRACTION_KEYS:
        out[key] = _spread([c[key] for c in cases if c.get(key) is not None])
    sxx = sum(x * x for x, _ in points)
    out["log_ratio_slope"] = _ratio(sum(x * y for x, y in points), sxx)
    out["n_allele_ratio_cases"] = len(points)
    out["span_offset_hist"] = hist
    out["c7_correct"] = {s: _ratio(ok, n) for s, (ok, n) in c7.items()}
    out["error_counts"] = _count_table(counts)
    out["c7_counts"] = c7
    return out
