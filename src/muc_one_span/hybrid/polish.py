"""S3/S7: POA draft, majority pileup polishing and homopolymer median vote.

Every tunable here (the POA sampling window, the insertion-vote majority fraction, the
homopolymer minimum run length, the number of rounds and the homopolymer-vote switch) is
a validated ``HybridSettings`` field (``muc_one_span.settings``). The parameters have no
defaults: callers pass the values of the ``RuntimeSettings`` they run with, so a user's
configuration can never be silently replaced by the package defaults.
"""

from __future__ import annotations

import random
import statistics
from collections import Counter
from typing import Any

from muc_one_span.hybrid.align import Columns, project
from muc_one_span.hybrid.poa import PoaBackend
from muc_one_span.hybrid.spans import SpanRead


def draft_consensus(
    members: list[SpanRead],
    n_poa: int,
    rng: random.Random,
    backend: PoaBackend,
    *,
    sample_window_floor_bp: float,
    sample_window_frac: float,
) -> str:
    """POA over a random sample of near-modal members (random, not quality-ranked).

    Members within ``max(sample_window_floor_bp, sample_window_frac * median)`` of the
    median length are "near"; the sample is drawn from those (falling back to all
    members when none are near, which cannot happen since the median member always
    qualifies).
    """
    median = statistics.median(m.length for m in members)
    tol = max(sample_window_floor_bp, sample_window_frac * median)
    near = [m for m in members if abs(m.length - median) <= tol] or list(members)
    sample = rng.sample(near, min(n_poa, len(near)))
    return backend.consensus([m.seq for m in sample])


def _projections(cons: str, full: list[str], partial: list[str]) -> list[tuple[str, Columns, bool]]:
    """Spanning reads align globally; fragments align as infixes and cover only part."""
    return [(r, project(r, cons), False) for r in full] + [
        (r, project(r, cons, partial=True), True) for r in partial
    ]


def pileup_polish(
    cons: str,
    full: list[str],
    partial: list[str] | None = None,
    *,
    insertion_majority_frac: float,
) -> tuple[str, int]:
    """One majority-vote round; each column/insertion slot is voted only by covering reads.

    A column takes its plurality base. An insertion slot is only accepted when its
    winning (nonempty) vote exceeds ``insertion_majority_frac`` of the covering votes.
    """
    n = len(cons)
    col: list[Counter[str | None]] = [Counter() for _ in range(n)]
    ins: list[Counter[str]] = [Counter() for _ in range(n + 1)]
    for _read, proj, is_partial in _projections(cons, full, partial or []):
        for pos in range(proj.first, proj.last):
            col[pos][proj.cols[pos]] += 1
        slots = range(proj.first + 1, proj.last) if is_partial else range(n + 1)
        for pos in slots:
            ins[pos][proj.ins.get(pos, "")] += 1
    out: list[str] = []
    changes = 0
    for pos in range(n + 1):
        if ins[pos]:
            ins_vote, ins_count = ins[pos].most_common(1)[0]
            if ins_vote and ins_count > insertion_majority_frac * sum(ins[pos].values()):
                out.append(ins_vote)
                changes += 1
        if pos < n:
            base = col[pos].most_common(1)[0][0] if col[pos] else cons[pos]
            if base != "-":
                out.append(str(base))
            changes += int(base != cons[pos])
    return "".join(out), changes


def _runs(seq: str, min_len: int) -> list[tuple[int, int, str]]:
    """Homopolymer runs of at least ``min_len`` as (start, end_excl, base)."""
    runs, i = [], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            runs.append((i, j, seq[i]))
        i = j
    return runs


def read_run_length(read: str, t2q: list[int], start: int, end: int, base: str) -> int:
    """Longest stretch of ``base`` in the read around the interval mapped to [start, end)."""
    qs, qe = t2q[start], t2q[end]
    a = qs
    while a > 0 and read[a - 1] == base:
        a -= 1
    b = max(qe, qs)
    while b < len(read) and read[b] == base:
        b += 1
    seg = "".join(c if c == base else " " for c in read[a:b])
    return max((len(m) for m in seg.split()), default=0)


def homopolymer_vote(
    cons: str,
    full: list[str],
    partial: list[str] | None = None,
    *,
    min_len: int,
) -> tuple[str, int]:
    """Set each run (>= min_len) to the median run length of the reads covering it."""
    runs = _runs(cons, min_len)
    if not runs or not (full or partial):
        return cons, 0
    projections = _projections(cons, full, partial or [])
    out, pos, changes = [], 0, 0
    for start, end, base in runs:
        lengths = [
            read_run_length(r, p.t2q, start, end, base)
            for r, p, _ in projections
            if p.covers(start, end)
        ]
        # +0.5 before truncating to int is the standard round-half-up rule (a format
        # definition, not a tunable); the result is floored at 1 base since a run cannot
        # have zero length once it exists.
        new_len = max(1, int(statistics.median(lengths) + 0.5)) if lengths else end - start
        out.append(cons[pos:start])
        out.append(base * new_len)
        changes += int(new_len != end - start)
        pos = end
    out.append(cons[pos:])
    return "".join(out), changes


def polish(
    cons: str,
    full: list[str],
    partial: list[str] | None = None,
    *,
    rounds: int,
    hp_vote: bool,
    insertion_majority_frac: float,
    hp_min_run: int,
) -> tuple[str, dict[str, Any]]:
    """Run ``rounds`` pileup rounds, each followed by an optional homopolymer vote."""
    info: dict[str, Any] = {"rounds": []}
    for _ in range(rounds):
        cons, changes = pileup_polish(
            cons, full, partial, insertion_majority_frac=insertion_majority_frac
        )
        hp = 0
        if hp_vote:
            cons, hp = homopolymer_vote(cons, full, partial, min_len=hp_min_run)
        info["rounds"].append({"changes": changes, "hp_changes": hp})
    return cons, info
