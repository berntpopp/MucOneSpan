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


def consensus_concordance(cons: str, full: list[str], partial: list[str] | None) -> float:
    """Mean per-base concordance between the final consensus and its covering reads.

    For every consensus position, this is the fraction of covering reads whose base
    there agrees with the consensus; the returned value is that fraction averaged over
    every position of ``cons``. This is real hybrid read-support evidence -- how well,
    on average, the polished consensus is backed by the reads that built/polished it --
    unlike the ladder's ``classify.py`` per-repeat ``confidence``, which only scores
    structural fit against the repeat dictionary, is fed the same way regardless of
    engine, and so says nothing about the hybrid engine's own reconstruction evidence
    (see ``hybrid/allele_fields.py::allele_info``, which reports this alongside a
    ``classification_confidence_status`` marker).

    An earlier version of this evidence required *unanimous* per-``unit_bp``-block
    agreement; at the read depths the hybrid engine routinely reconstructs from
    (dozens to thousands of spanning reads), realistic per-base sequencing noise means
    at least one covering read disagrees somewhere in nearly every block, collapsing
    that stricter definition to ~0 regardless of how good the consensus actually is
    (confirmed on real reconstructed alleles). A per-position mean does not have that
    problem: it stays informative at any depth. A position with no covering read
    contributes 0 (fail-closed), matching the rest of this project's evidence handling.
    An empty consensus is that same "no evidence" case in the limit -- there is no
    reconstruction for any read to support -- so it also fails closed to 0.0 rather
    than a misleadingly perfect 1.0.
    """
    n = len(cons)
    if n == 0:
        return 0.0
    match = [0] * n
    cover = [0] * n
    for _read, proj, _is_partial in _projections(cons, full, partial or []):
        for pos in range(proj.first, proj.last):
            cover[pos] += 1
            if proj.cols[pos] == cons[pos]:
                match[pos] += 1
    return sum(match[i] / cover[i] if cover[i] else 0.0 for i in range(n)) / n


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
