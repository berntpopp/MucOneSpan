"""S1: locate motif-1/motif-9 anchors, orient reads, and categorise them."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import infix_hit, rc
from muc_one_span.settings import HybridSettings

FLANK_ANCHOR_BP = 30
UNIT = 60


@dataclass(frozen=True)
class ReadRecord:
    """One input read as stored in FASTQ."""

    name: str
    seq: str
    qual: str


@dataclass(frozen=True)
class SpanRead:
    """A spanning read trimmed to motif 1..motif 9 and oriented to the VNTR forward strand."""

    name: str
    seq: str
    mean_q: float
    strand: str  # "+" read was already forward, "-" read was reverse-complemented
    anchor_edits: int
    anchor_basis: str  # "motif" or "flank"

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass(frozen=True)
class Anchors:
    """Motif 1/9 anchors plus flank anchors used when a motif carries a mutation."""

    left: str
    right: str
    left_flank: str
    right_flank: str

    @classmethod
    def from_dictionary(cls, rd: RepeatDictionary) -> Anchors:
        return cls(
            rd.repeats["1"],
            rd.repeats["9"],
            rd.flanking_left[-FLANK_ANCHOR_BP:],
            rd.flanking_right[:FLANK_ANCHOR_BP],
        )


@dataclass
class ReadCategories:
    """Reads grouped by anchor evidence."""

    spanning: list[SpanRead] = field(default_factory=list)
    left_anchored: list[ReadRecord] = field(default_factory=list)
    right_anchored: list[ReadRecord] = field(default_factory=list)
    internal_or_offtarget: list[ReadRecord] = field(default_factory=list)

    def counts(self) -> dict[str, int]:
        return {
            "spanning": len(self.spanning),
            "left_anchored": len(self.left_anchored),
            "right_anchored": len(self.right_anchored),
            "internal_or_offtarget": len(self.internal_or_offtarget),
        }


def _mean_q(qual: str) -> float:
    return sum(ord(c) - 33 for c in qual) / len(qual) if qual else 0.0


def _span_in(target: str, anchors: Anchors, k: int) -> tuple[int, int, int, str] | None:
    """Return (start, end_excl, edits, basis) of motif1..motif9 in target, or None."""
    left = infix_hit(anchors.left, target, k)
    right = None
    if left is not None:
        tail = infix_hit(anchors.right, target[left[1] :], k)
        right = None if tail is None else (tail[0] + left[1], tail[1] + left[1], tail[2])
    if left is not None and right is not None:
        return left[0], right[1], left[2] + right[2], "motif"
    kf = max(2, k // 4)
    lf = infix_hit(anchors.left_flank, target, kf)
    rf = infix_hit(anchors.right_flank, target[lf[1] :], kf) if lf else None
    if lf is None or rf is None:
        return None
    return lf[1], rf[0] + lf[1], lf[2] + rf[2], "flank"


def categorize_reads(
    reads: Iterable[ReadRecord], anchors: Anchors, settings: HybridSettings
) -> ReadCategories:
    """Classify reads as spanning (oriented, trimmed) / one-end anchored / other."""
    k = settings.anchor_max_edits
    lo, hi = settings.min_span_units * UNIT, settings.max_span_units * UNIT
    cats = ReadCategories()
    for rec in reads:
        seq = rec.seq.upper()
        best: tuple[tuple[int, int, int, str], str, str, str] | None = None
        for strand, target, qual in (("+", seq, rec.qual), ("-", rc(seq), rec.qual[::-1])):
            hit = _span_in(target, anchors, k)
            if hit and (best is None or hit[2] < best[0][2]):
                best = (hit, strand, target, qual)
        if best is not None:
            (start, end, edits, basis), strand, target, qual = best
            if lo <= end - start <= hi:
                cats.spanning.append(
                    SpanRead(
                        rec.name, target[start:end], _mean_q(qual[start:end]), strand, edits, basis
                    )
                )
                continue
        has_left = any(infix_hit(anchors.left, t, k) for t in (seq, rc(seq)))
        has_right = any(infix_hit(anchors.right, t, k) for t in (seq, rc(seq)))
        if has_left and not has_right:
            cats.left_anchored.append(rec)
        elif has_right and not has_left:
            cats.right_anchored.append(rec)
        else:
            cats.internal_or_offtarget.append(rec)
    return cats
