"""S8/S10: residual heterogeneity QC and per-event read-level support.

Producer contract for ``read_support.status == "supported"`` (consumed unchanged by
``clinical_gates.mutation_supported``); every threshold is a ``HybridSettings`` field:

* Homopolymer event (dictionary template = single-base indel inside a consensus run
  ``>= hp_event_min_run``): ``n >= hp_min_reads``, stutter-aware LLR ``>= hp_llr_min``,
  alt fraction ``>= hp_min_alt_frac``, and no strand with ``>= hp_min_strand_reads``
  reads has a negative strand LLR (else ``discordant``). A strand with zero reads
  never fails the event.
* Other templated events: parent-vs-template edit-distance competition per read;
  ``n >= hp_min_reads``, ``alt / n >= hp_min_alt_frac`` and ``alt > ref``.
* ``n < hp_min_reads`` is ``insufficient_depth``: blocked, never negative.
* A mutation whose repeat unit or parent type cannot be located is ``not_localized``.

Only values in ``clinical_gates.READ_SUPPORT_STATUSES`` are emitted.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import Columns, edit_distance, global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.settings import HybridSettings

Profile = dict[str, list[float]]
STRANDS = ("+", "-")
# Display precision of reported fractions and log-likelihood ratios (output format
# only: every status is decided on the unrounded values).
FRACTION_DECIMALS = 3
LLR_DECIMALS = 1


def residual_sites(cons: str, reads: list[str], af: float, *, min_run: int) -> list[dict[str, Any]]:
    """Consensus columns where one non-consensus base/deletion reaches ``af``.

    Columns inside consensus homopolymer runs ``>= min_run`` are skipped: an indel in
    a run has no unique column, so its per-column frequency is not a mixture signal.
    """
    in_run = {i for s, e, _ in _runs(cons, min_run) for i in range(s, e)}
    counts: list[Counter[str | None]] = [Counter() for _ in cons]
    for read in reads:
        for pos, base in enumerate(global_columns(read, cons).cols):
            counts[pos][base] += 1
    out = []
    for pos, c in enumerate(counts):
        if pos in in_run or not c:
            continue
        tot = sum(c.values())
        for base, n in c.most_common():
            if base != cons[pos] and n / tot >= af:
                af_obs = round(n / tot, FRACTION_DECIMALS)
                out.append({"pos": pos, "ref": cons[pos], "alt": base, "af": af_obs, "n": tot})
                break
    return out


def _smooth(counter: Counter[int], s: HybridSettings) -> list[float]:
    size = s.hp_max_run_len + 1
    tot = sum(counter.values())
    pseudo = s.hp_background_pseudocount
    return [(counter.get(i, 0) + pseudo) / (tot + pseudo * size) for i in range(size)]


def _shift(p: list[float], d: int) -> list[float]:
    cap = len(p) - 1
    q = [p[min(max(i - d, 0), cap)] for i in range(len(p))]
    total = sum(q)
    return [x / total for x in q]


def homopolymer_background(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    base: str,
    length: int,
    exclude: int,
    s: HybridSettings,
) -> Profile:
    """Per-strand observed-length profile at consensus runs of ``base`` x ``length``.

    The event run (the one containing ``exclude``) is left out; each profile has
    ``hp_max_run_len + 1`` smoothed entries (longer observations are capped).
    """
    runs = [
        (st, e)
        for st, e, b in _runs(cons, length)
        if b == base and e - st == length and not st <= exclude < e
    ]
    per: dict[str, Counter[int]] = {strand: Counter() for strand in STRANDS}
    for seq, strand, proj in reads:
        for st, e in runs:
            observed = read_run_length(seq, proj.t2q, st, e, base)
            per.setdefault(strand, Counter())[min(observed, s.hp_max_run_len)] += 1
    return {strand: _smooth(c, s) for strand, c in per.items()}


def homopolymer_llr(obs: list[tuple[str, int]], background: Profile, shift: int) -> float:
    """log L(event length) - log L(reference length) over (strand, observed length).

    The reference model is the strand's background profile; the event model is that
    profile shifted by ``shift`` bases (the template's run-length change, +1 for a
    single-base insertion and -1 for a deletion). An observation on a strand without a
    background profile carries no information and contributes 0.
    """
    llr = 0.0
    for strand, observed in obs:
        p0 = background.get(strand)
        if not p0:
            continue
        p1 = _shift(p0, shift)
        k = min(observed, len(p0) - 1)
        llr += math.log(p1[k] / p0[k])
    return llr


def homopolymer_event_run(
    mutation: dict[str, Any],
    rd: RepeatDictionary,
    cons: str,
    start: int,
    end: int,
    s: HybridSettings,
) -> tuple[int, int, str, int] | None:
    """(run_start, run_end, base, shift) when the template is a homopolymer event.

    The type comes from the dictionary template (C8.1), never from the consensus: a
    single insert of one base, or a single one-base delete, that lies in a consensus
    run ``>= hp_event_min_run`` inside the unit ``[start, end)``. Runs whose event or
    reference length reaches ``hp_max_run_len`` are not modelled (returns None, so the
    event falls back to parent-vs-template competition).
    """
    template = rd.mutations.get(mutation.get("mutation_name") or "") or {}
    changes = template.get("changes") or []
    if len(changes) != 1:
        return None
    change = changes[0]
    offset = int(change["start"]) - 1  # dictionary coordinates are 1-based
    if change.get("type") == "insert" and len(change.get("sequence", "")) == 1:
        base, shift = str(change["sequence"]).upper(), 1
    elif change.get("type") == "delete" and int(change["end"]) == int(change["start"]):
        parent = rd.repeats.get(mutation.get("closest_type") or "")
        if parent is None or not 0 <= offset < len(parent):
            return None
        base, shift = parent[offset].upper(), -1
    else:
        return None
    pos = start + offset
    for rs, re_, b in _runs(cons[start:end], s.hp_event_min_run):
        run_start, run_end = rs + start, re_ + start
        if b == base and run_start <= pos <= run_end:
            if max(run_end - run_start, run_end - run_start - shift) >= s.hp_max_run_len:
                return None
            return run_start, run_end, b, shift
    return None


def hp_status(
    n: int,
    llr: float,
    alt_frac: float,
    strand_llr: dict[str, float],
    strand_n: dict[str, int],
    s: HybridSettings,
) -> str:
    """Status of a homopolymer event (spec §5 thresholds plus strand consistency)."""
    if n < s.hp_min_reads:
        return "insufficient_depth"
    if llr < s.hp_llr_min or alt_frac < s.hp_min_alt_frac:
        return "not_supported"
    for strand, count in strand_n.items():
        if count > 0 and count >= s.hp_min_strand_reads and strand_llr.get(strand, 0.0) < 0:
            return "discordant"
    return "supported"


def competition_status(n: int, alt: int, ref: int, s: HybridSettings) -> str:
    """Status of a parent-vs-template competition event (C8.3)."""
    if n < s.hp_min_reads:
        return "insufficient_depth"
    if alt / n >= s.hp_min_alt_frac and alt > ref:
        return "supported"
    return "not_supported"


def _strand_fracs(flags: dict[str, list[bool]]) -> dict[str, float | None]:
    return {
        strand: round(sum(v) / len(v), FRACTION_DECIMALS) if v else None
        for strand, v in flags.items()
    }


def _homopolymer_support(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    run: tuple[int, int, str, int],
    s: HybridSettings,
) -> dict[str, Any]:
    start, end, base, shift = run
    length = end - start
    background = homopolymer_background(cons, reads, base, length - shift, start, s)
    obs = [(st, read_run_length(seq, p.t2q, start, end, base)) for seq, st, p in reads]

    def is_alt(k: int) -> bool:
        return k >= length if shift > 0 else k <= length

    n = len(obs)
    alt = sum(1 for _, k in obs if is_alt(k))
    strands = sorted({*STRANDS, *(st for st, _ in obs)})
    strand_obs = {st: [o for o in obs if o[0] == st] for st in strands}
    strand_llr = {st: homopolymer_llr(v, background, shift) for st, v in strand_obs.items()}
    strand_n = {st: len(v) for st, v in strand_obs.items()}
    llr = homopolymer_llr(obs, background, shift)
    alt_frac = alt / n if n else 0.0
    return {
        "kind": "homopolymer",
        "n": n,
        "alt": alt,
        "ref": n - alt,
        "other": 0,
        "alt_frac": round(alt_frac, FRACTION_DECIMALS),
        "strand_alt_frac": _strand_fracs(
            {st: [is_alt(k) for _, k in v] for st, v in strand_obs.items()}
        ),
        "llr": round(llr, LLR_DECIMALS),
        "strand_llr": {k: round(v, LLR_DECIMALS) for k, v in strand_llr.items()},
        "status": hp_status(n, llr, alt_frac, strand_llr, strand_n, s),
    }


def _competition_support(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    start: int,
    end: int,
    parent: str,
    s: HybridSettings,
) -> dict[str, Any]:
    window = cons[start:end]
    alt = ref = other = 0
    per_strand: dict[str, list[bool]] = {strand: [] for strand in STRANDS}
    for seq, strand, proj in reads:
        piece = seq[proj.t2q[start] : proj.t2q[end]]
        d_alt, d_ref = edit_distance(piece, window), edit_distance(piece, parent)
        alt += d_alt < d_ref
        ref += d_alt > d_ref
        other += d_alt == d_ref
        per_strand.setdefault(strand, []).append(d_alt < d_ref)
    n = alt + ref + other
    return {
        "kind": "competition",
        "n": n,
        "alt": alt,
        "ref": ref,
        "other": other,
        "alt_frac": round(alt / n if n else 0.0, FRACTION_DECIMALS),
        "strand_alt_frac": _strand_fracs(per_strand),
        "status": competition_status(n, alt, ref, s),
    }


def _not_localized() -> dict[str, Any]:
    return {
        "kind": "none",
        "n": 0,
        "alt": 0,
        "ref": 0,
        "other": 0,
        "alt_frac": 0.0,
        "strand_alt_frac": dict.fromkeys(STRANDS),
        "status": "not_localized",
    }


def event_read_support(
    cons: str,
    classification: dict[str, Any],
    reads: list[tuple[str, str]],
    rd: RepeatDictionary,
    settings: HybridSettings,
) -> dict[int, dict[str, Any]]:
    """Read-level support per detected mutation (keyed by its index in the classification).

    ``reads`` are (sequence, strand) pairs of spanning reads assigned to this allele,
    already oriented to the consensus.
    """
    repeats = {r.get("index"): r for r in classification.get("repeats", [])}
    projected = [(seq, strand, global_columns(seq, cons)) for seq, strand in reads]
    results: dict[int, dict[str, Any]] = {}
    for idx, mut in enumerate(classification.get("mutations_detected", [])):
        rep = repeats.get(mut.get("repeat_index"))
        parent = rd.repeats.get(mut.get("closest_type") or "")
        if rep is None or parent is None:
            results[idx] = _not_localized()
            continue
        start, end = int(rep["start"]), int(rep["end"])
        run = homopolymer_event_run(mut, rd, cons, start, end, settings)
        results[idx] = (
            _homopolymer_support(cons, projected, run, settings)
            if run is not None
            else _competition_support(cons, projected, start, end, parent, settings)
        )
    return results
