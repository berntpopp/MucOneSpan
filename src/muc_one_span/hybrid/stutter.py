"""Run-length-aware homopolymer stutter background (Task 15f).

Stutter depends on the run length: a longer run is misread by one base more often.
A homopolymer event is scored as a mixture of two alleles, the event run length and
the no-event run length (``evidence._homopolymer_support``), and each allele is
convolved with the per-strand stutter profile of **its own** length.

A profile is an observed-length distribution (index = observed run length, capped at
``hp_max_run_len``), smoothed with ``hp_background_pseudocount``. For one base and one
strand, the profile of a length is chosen in this order (``hp_stutter_model =
"length"``):

1. **Measured**: the observations at the sample's own peer runs of that base and
   length (the event run is left out), when there are at least
   ``hp_stutter_min_class_runs`` such runs and ``hp_stutter_min_class_reads`` clean
   observations on the strand.
2. **Extrapolated** from the two measured lengths nearest to it (ties to the shorter):
   each observed-minus-true error value follows its own geometric trend in run
   length (log-linear), with a per-base growth factor capped to
   ``[1 / hp_stutter_max_growth, hp_stutter_max_growth]``, and the result is
   renormalised. Between two measured lengths this interpolates.
3. **Shifted**: the single measured length's profile moved to this length.
4. Otherwise the pre-15f rule: the no-event allele uses whatever its own length's
   peer runs show, and the event allele that profile shifted by the event.

Lengths are per base; there is no pooling across bases. ``hp_stutter_model =
"shift"`` is the pre-15f model (rule 4 for both alleles). Every tunable is a
validated ``HybridSettings`` field.
"""

from __future__ import annotations

from collections import Counter

from muc_one_span.hybrid.align import Columns
from muc_one_span.hybrid.polish import _runs, run_observation
from muc_one_span.settings import HybridSettings

Profile = dict[str, list[float]]
STRANDS = ("+", "-")
# A homopolymer run has at least two bases (definition, not a tunable).
MIN_RUN = 2


def smooth(counter: Counter[int], s: HybridSettings) -> list[float]:
    """Observed-length distribution with ``hp_background_pseudocount`` per length."""
    size = s.hp_max_run_len + 1
    tot = sum(counter.values())
    pseudo = s.hp_background_pseudocount
    return [(counter.get(i, 0) + pseudo) / (tot + pseudo * size) for i in range(size)]


def shift(p: list[float], d: int) -> list[float]:
    """Move a profile by ``d`` bases (edge entries repeat), renormalised."""
    cap = len(p) - 1
    q = [p[min(max(i - d, 0), cap)] for i in range(len(p))]
    total = sum(q)
    return [x / total for x in q]


def class_counts(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    base: str,
    length: int,
    exclude: int,
    s: HybridSettings,
) -> tuple[int, dict[str, Counter[int]]]:
    """(peer runs, per-strand observed lengths) at consensus runs of ``base`` x ``length``.

    The run containing ``exclude`` (the event run) is left out, and so is every read
    that does not observe a run cleanly (``run_observation``); observations above
    ``hp_max_run_len`` are capped.
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
    return len(runs), per


def length_profile(
    measured: dict[int, list[float]], length: int, max_growth: float
) -> list[float] | None:
    """Profile of ``length`` from the measured lengths (rules 1-3), None without any.

    Extrapolation aligns the two nearest measured profiles by error value (each moved
    to ``length``) and grows each error value's share geometrically from the second to
    the first, at a per-base factor capped to ``[1 / max_growth, max_growth]``. A share
    that is zero at the nearest length stays zero; one that is zero only at the second
    grows at the cap.
    """
    if length in measured:
        return measured[length]
    near = sorted(measured, key=lambda m: (abs(m - length), m))
    if not near:
        return None
    first = shift(measured[near[0]], length - near[0])
    if len(near) == 1:
        return first
    second = shift(measured[near[1]], length - near[1])
    span = near[0] - near[1]
    raw = []
    for a, b in zip(first, second, strict=True):
        if a <= 0:
            raw.append(0.0)
            continue
        # a / b is infinite for b = 0: the cap in the direction of the trend.
        per_base = (a / b) ** (1 / span) if b > 0 else max_growth if span > 0 else 0.0
        per_base = min(max(per_base, 1 / max_growth), max_growth)
        raw.append(a * per_base ** (length - near[0]))
    total = sum(raw)
    return [x / total for x in raw]


def allele_profiles(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    run: tuple[int, int, str, int],
    strands: set[str],
    s: HybridSettings,
) -> tuple[Profile, Profile]:
    """(event, no-event) per-strand profiles of a homopolymer event run.

    ``run`` is ``(start, end, base, shift)`` from ``evidence.homopolymer_event_run``;
    the event allele has the consensus run length and the no-event allele that length
    minus ``shift``. Profiles are built for ``STRANDS`` and every strand in
    ``strands``; only strands in ``strands`` (those with event observations) search
    for measured lengths.
    """
    start, end, base, d = run
    event_len = end - start
    none_len = event_len - d
    cache: dict[int, tuple[int, dict[str, Counter[int]]]] = {}

    def counts(length: int) -> tuple[int, dict[str, Counter[int]]]:
        if length not in cache:
            cache[length] = class_counts(cons, reads, base, length, start, s)
        return cache[length]

    def measured(length: int, strand: str) -> list[float] | None:
        n_runs, per = counts(length)
        c = per.get(strand, Counter())
        if n_runs < s.hp_stutter_min_class_runs or sum(c.values()) < s.hp_stutter_min_class_reads:
            return None
        return smooth(c, s)

    def profile(strand: str, length: int) -> list[float] | None:
        found: dict[int, list[float]] = {}
        order = sorted(range(MIN_RUN, s.hp_max_run_len + 1), key=lambda m: (abs(m - length), m))
        for m in order:
            p = measured(m, strand)
            if p is not None:
                found[m] = p
                if m == length or len(found) == 2:
                    break
        return length_profile(found, length, s.hp_stutter_max_growth)

    _n, per_none = counts(none_len)
    names = sorted({*STRANDS, *per_none, *strands})
    none = {st: smooth(per_none.get(st, Counter()), s) for st in names}
    event = {st: shift(p, d) for st, p in none.items()}
    if s.hp_stutter_model == "shift":
        return event, none
    for st in sorted(strands):
        none[st] = profile(st, none_len) or none[st]
        event[st] = profile(st, event_len) or shift(none[st], d)
    return event, none
