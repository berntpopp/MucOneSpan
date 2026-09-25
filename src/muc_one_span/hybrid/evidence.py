"""S8/S10: residual heterogeneity QC and per-event read-level support.

Producer contract for ``read_support.status == "supported"`` (consumed unchanged by
``clinical_gates.mutation_supported``); every threshold is a ``HybridSettings`` field:

* Homopolymer event (dictionary template = single-base indel inside a consensus run
  ``>= hp_event_min_run``): ``n >= hp_min_reads``, stutter-aware LLR ``>= hp_llr_min``,
  alt fraction ``>= hp_min_alt_frac``, the no-event share of the maximum-likelihood
  event/no-event stutter mixture ``<= event_max_alternative_frac`` (else
  ``discordant``), and no strand with ``>= hp_min_strand_reads`` reads has a negative
  strand LLR (else ``discordant``). A strand with zero reads never fails the event.
  Only reads that keep both consensus bases bounding the run, with nothing but the
  run's base between them, observe its length; any other read is ``other`` (never
  support).
* Other events: per-read edit-distance competition over the event unit extended by
  ``event_context_units`` repeat units on each side, between the event allele (the
  consensus) and the best alternative, which is the no-event allele (the unit reverted
  to its dictionary parent) or the read-derived allele (the polish-rule pileup
  consensus of the reads that do not favour the event). A tie is ``other``. ``n >= hp_min_reads``,
  ``alt / n >= hp_min_alt_frac`` and ``alt > ref``, and the alternative share
  ``ref / n <= event_max_alternative_frac`` (else ``discordant``).
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
from muc_one_span.hybrid.polish import _runs, polish, run_observation
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

    The event run (the one containing ``exclude``) is left out, and so is every read
    that does not observe a run cleanly (``run_observation``); each profile has
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
            observed = run_observation(seq, proj, cons, st, e)
            if observed is not None:
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
    # Runs are found in the whole consensus so a run that crosses the unit boundary is
    # measured whole; only runs that intersect the unit and contain the event count.
    for run_start, run_end, b in _runs(cons, s.hp_event_min_run):
        if run_end <= start or run_start >= end:
            continue
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
    *,
    alternative_frac: float,
) -> str:
    """Status of a homopolymer event (spec §5 thresholds plus strand consistency).

    ``alternative_frac`` is the no-event share of the event/no-event mixture
    (``1 - event_allele_fraction``); above ``event_max_alternative_frac`` the reads are a
    mixture rather than the event allele, so the event is ``discordant``. Zero reads is
    always ``insufficient_depth``, whatever the configured minimum.
    """
    if n <= 0 or n < s.hp_min_reads:
        return "insufficient_depth"
    if llr < s.hp_llr_min or alt_frac < s.hp_min_alt_frac:
        return "not_supported"
    if alternative_frac > s.event_max_alternative_frac:
        return "discordant"
    for strand, count in strand_n.items():
        if count > 0 and count >= s.hp_min_strand_reads and strand_llr.get(strand, 0.0) < 0:
            return "discordant"
    return "supported"


def event_allele_fraction(obs: list[tuple[float, float]]) -> float:
    """Maximum-likelihood weight of the event allele in an event/no-event read mixture.

    Each observation is ``(p_event, p_no_event)``: the probability of one read's
    observation under each allele. The weight ``f`` maximises
    ``sum(log(f * p_event + (1 - f) * p_no_event))`` on ``[0, 1]``. That log-likelihood
    is concave, so its slope decreases in ``f`` and bisection on the slope's sign finds
    the maximum; the loop ends when the float interval cannot shrink any further (no
    tolerance tunable). No informative observation returns 0.0 (fail closed).
    """

    def slope(f: float) -> float:
        total = 0.0
        for p1, p0 in obs:
            den = f * p1 + (1 - f) * p0
            if den <= 0:  # a read impossible at f (e.g. a column allele at f = 0 or 1)
                return math.inf if p1 > p0 else -math.inf
            total += (p1 - p0) / den
        return total

    if not obs or slope(0.0) <= 0:
        return 0.0
    if slope(1.0) >= 0:
        return 1.0
    lo, hi = 0.0, 1.0
    while True:
        mid = (lo + hi) / 2
        if mid <= lo or mid >= hi:
            return mid
        if slope(mid) > 0:
            lo = mid
        else:
            hi = mid


def competition_status(n: int, alt: int, ref: int, s: HybridSettings) -> str:
    """Status of an event-vs-alternative competition event (C8.3); zero reads is blocked.

    ``ref`` counts reads that favour an alternative over the event allele; a share
    ``ref / n`` above ``event_max_alternative_frac`` is ``discordant``.
    """
    if n <= 0 or n < s.hp_min_reads:
        return "insufficient_depth"
    if alt / n < s.hp_min_alt_frac or alt <= ref:
        return "not_supported"
    if ref / n > s.event_max_alternative_frac:
        return "discordant"
    return "supported"


def _strand_fracs(flags: dict[str, list[bool]]) -> dict[str, float | None]:
    return {
        strand: round(sum(v) / len(v), FRACTION_DECIMALS) if v else None
        for strand, v in flags.items()
    }


def _mixture_obs(
    obs: list[tuple[str, int]], background: Profile, shift: int
) -> list[tuple[float, float]]:
    """(p_event, p_no_event) per observation; a strand without a profile is skipped."""
    out = []
    for strand, observed in obs:
        p0 = background.get(strand)
        if not p0:
            continue
        k = min(observed, len(p0) - 1)
        out.append((_shift(p0, shift)[k], p0[k]))
    return out


def _homopolymer_support(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    run: tuple[int, int, str, int],
    s: HybridSettings,
) -> dict[str, Any]:
    start, end, base, shift = run
    length = end - start
    background = homopolymer_background(cons, reads, base, length - shift, start, s)
    observed = [(st, run_observation(seq, p, cons, start, end)) for seq, st, p in reads]
    obs = [(st, k) for st, k in observed if k is not None]

    def is_alt(k: int) -> bool:
        return k >= length if shift > 0 else k <= length

    n = len(reads)
    alt = sum(1 for _, k in obs if is_alt(k))
    strands = sorted({*STRANDS, *(st for st, _ in obs)})
    strand_obs = {st: [o for o in obs if o[0] == st] for st in strands}
    strand_llr = {st: homopolymer_llr(v, background, shift) for st, v in strand_obs.items()}
    strand_n = {st: len(v) for st, v in strand_obs.items()}
    llr = homopolymer_llr(obs, background, shift)
    alt_frac = alt / n if n else 0.0
    alternative_frac = 1.0 - event_allele_fraction(_mixture_obs(obs, background, shift))
    return {
        "kind": "homopolymer",
        "n": n,
        "alt": alt,
        "ref": len(obs) - alt,
        "other": n - len(obs),
        "alt_frac": round(alt_frac, FRACTION_DECIMALS),
        "alternative_frac": round(alternative_frac, FRACTION_DECIMALS),
        "strand_alt_frac": _strand_fracs(
            {st: [is_alt(k) for _, k in v] for st, v in strand_obs.items()}
        ),
        "llr": round(llr, LLR_DECIMALS),
        "strand_llr": {k: round(v, LLR_DECIMALS) for k, v in strand_llr.items()},
        "status": hp_status(
            n, llr, alt_frac, strand_llr, strand_n, s, alternative_frac=alternative_frac
        ),
    }


def _favours(piece: str, event: str, alternatives: list[str]) -> int:
    """+1 when the read piece is closer to the event allele than to every alternative,
    -1 when an alternative is closer, 0 for a tie (a tie is never support)."""
    d_event = edit_distance(piece, event)
    d_alternative = min(edit_distance(piece, a) for a in alternatives)
    return (d_event < d_alternative) - (d_alternative < d_event)


def _pileup_consensus(template: str, pieces: list[str], s: HybridSettings) -> str:
    return polish(
        template,
        pieces,
        rounds=s.polish_rounds,
        hp_vote=s.hp_vote,
        insertion_majority_frac=s.polish_insertion_majority_frac,
        hp_min_run=s.hp_vote_min_run,
    )[0]


def _competition_support(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    start: int,
    end: int,
    parent: str,
    context_bp: int,
    s: HybridSettings,
) -> dict[str, Any]:
    lo, hi = max(0, start - context_bp), min(len(cons), end + context_bp)
    event = cons[lo:hi]
    no_event = cons[lo:start] + parent + cons[end:hi]
    pieces = [(seq[proj.t2q[lo] : proj.t2q[hi]], strand) for seq, strand, proj in reads]
    # Best read-derived alternative: the pileup/homopolymer-vote consensus (the polish
    # stage's own rules and settings) of the reads that do not favour the event over
    # the no-event allele. It is the allele those reads carry even when the consensus
    # error sits outside the unit the classifier blamed (a displaced indel).
    first = [_favours(piece, event, [no_event]) for piece, _ in pieces]
    rest = [piece for (piece, _), vote in zip(pieces, first, strict=True) if vote <= 0 and piece]
    derived = _pileup_consensus(event, rest, s) if rest else event
    alternatives = [no_event] + ([derived] if derived not in (event, no_event) else [])
    alt = ref = other = 0
    per_strand: dict[str, list[bool]] = {strand: [] for strand in STRANDS}
    for piece, strand in pieces:
        vote = _favours(piece, event, alternatives)
        alt += vote > 0
        ref += vote < 0
        other += vote == 0
        per_strand.setdefault(strand, []).append(vote > 0)
    n = alt + ref + other
    return {
        "kind": "competition",
        "n": n,
        "alt": alt,
        "ref": ref,
        "other": other,
        "alt_frac": round(alt / n if n else 0.0, FRACTION_DECIMALS),
        "alternative_frac": round(ref / n if n else 0.0, FRACTION_DECIMALS),
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
        "alternative_frac": None,
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
    # Context around each unit, in repeat units of the dictionary's unit length.
    context_bp = round(settings.event_context_units * rd.repeat_length_bp)
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
            else _competition_support(cons, projected, start, end, parent, context_bp, settings)
        )
    return results
