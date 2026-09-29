"""S2 PCR dimer products: spanning reads that hold two amplicons joined head to tail.

A head-to-tail PCR/ligation dimer of amplicons a and b is read as
``flank + allele_a + flank + flank + allele_b + flank``. The S1 span search takes the
best motif-1 and the best motif-9 hit after it, so when they fall in different copies
the "span" is ``allele_a + amplicon junction + allele_b``: a length peak near
``L_a + L_b`` (``2 x L`` for a homozygous dimer) that no allele explains.

A read is a dimer product only on structural evidence, never on length alone: inside
the span (the terminal motif-1 and motif-9 units excluded) it carries a motif-9 hit
followed by a motif-1 hit (the amplicon junction), and the two parts it splits into
each fall in an accepted allele peak's assignment window. Motif 1 is far from every
other repeat unit of the dictionary, so a genuine allele read, whatever its length, has
no such junction: a real allele at ``2 x L`` is never taken for a dimer.

Dimer reads of one parent pair are accepted as an artefact only while they stay a
minority (``dimer_max_parent_frac`` of the smaller parent's support). Otherwise the pair
is left unexplained and its reads are evaluated like any other length candidate.
"""

from __future__ import annotations

import statistics
from collections.abc import Sequence
from typing import Any, NamedTuple

from muc_one_span.hybrid.align import infix_hit
from muc_one_span.hybrid.spans import Anchors, SpanRead
from muc_one_span.settings import HybridSettings


class Junction(NamedTuple):
    """Where a dimer read splits: part lengths (motif 1..9 each) and the junction gap."""

    first_bp: int
    second_bp: int
    gap_bp: int


class DimerCall(NamedTuple):
    """Reads explained as dimer products, one rejected-peak record per parent pair, and
    the centres of the accepted peaks that are parents of an explained pair."""

    reads: list[SpanRead]
    entries: list[dict[str, Any]]
    parents: list[float]


def find_junction(span: SpanRead, anchors: Anchors, settings: HybridSettings) -> Junction | None:
    """The internal motif-9 -> motif-1 junction of a spanning read, or None.

    The search skips the span's own terminal units (its motif 1 and motif 9) and uses
    the S1 anchor edit allowance (``anchor_max_edits``).
    """
    unit = anchors.unit_bp
    inner = span.seq[unit:-unit]
    k = settings.anchor_max_edits
    end_of_first = infix_hit(anchors.right, inner, k)
    if end_of_first is None:
        return None
    start_of_second = infix_hit(anchors.left, inner[end_of_first[1] :], k)
    if start_of_second is None:
        return None
    first = unit + end_of_first[1]
    second_start = first + start_of_second[0]
    return Junction(first, span.length - second_start, start_of_second[0])


def _parent(length: int, peaks: Sequence[tuple[float, float]]) -> int | None:
    """Index of the accepted peak whose assignment window holds ``length``, if any."""
    for index, (center, half_window) in enumerate(peaks):
        if abs(length - center) <= half_window:
            return index
    return None


def recognise_dimers(
    candidates: Sequence[SpanRead],
    peaks: Sequence[tuple[float, float, int]],
    anchors: Anchors,
    settings: HybridSettings,
) -> DimerCall:
    """Pick the dimer products out of ``candidates``.

    ``peaks`` are the accepted peaks as ``(center_bp, half_window_bp, support)``. A read
    is a dimer product of the pair ``(i, j)`` when its junction parts fall in the
    windows of peaks ``i`` and ``j``; a pair is accepted when its read count is at most
    ``dimer_max_parent_frac`` times the smaller parent's support.
    """
    if not peaks:
        return DimerCall([], [], [])
    windows = [(center, half) for center, half, _ in peaks]
    shortest = min(center - half for center, half in windows)
    by_pair: dict[tuple[int, int], list[SpanRead]] = {}
    for span in candidates:
        if span.length < 2 * shortest:
            continue  # shorter than any two accepted alleles: cannot be a dimer
        junction = find_junction(span, anchors, settings)
        if junction is None:
            continue
        i, j = _parent(junction.first_bp, windows), _parent(junction.second_bp, windows)
        if i is None or j is None:
            continue
        by_pair.setdefault((min(i, j), max(i, j)), []).append(span)
    reads: list[SpanRead] = []
    entries: list[dict[str, Any]] = []
    parents: set[float] = set()
    unit = anchors.unit_bp
    for (i, j), members in sorted(by_pair.items()):
        parent_support = min(peaks[i][2], peaks[j][2])
        if len(members) > settings.dimer_max_parent_frac * parent_support:
            continue  # not a minority artefact: leave it unexplained
        reads.extend(members)
        parents |= {peaks[i][0], peaks[j][0]}
        center = statistics.median(s.length for s in members)
        entries.append(
            {
                "center_bp": round(center, 1),
                "units": round(center / unit),
                "support": len(members),
                "reason": "dimer",
                "parent_units": [round(peaks[i][0] / unit), round(peaks[j][0] / unit)],
            }
        )
    return DimerCall(reads, entries, sorted(parents))
