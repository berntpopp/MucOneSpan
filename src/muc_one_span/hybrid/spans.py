"""S1: locate motif-1/motif-9 anchors, orient reads, and categorise them."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import infix_hit, rc
from muc_one_span.settings import HybridSettings, ReferenceLayoutSettings

PHRED_OFFSET = 33  # FASTQ Phred+33 quality encoding (format definition, not a tunable)


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
    """Motif 1/9 anchors plus flank anchors used when a motif carries a mutation.

    ``unit_bp`` is the repeat-unit length taken from the loaded dictionary
    (``rd.repeat_length_bp``, 60 bp for the bundled MUC1 dictionary); it is the
    single source of the "unit" used to scale spans and length-model windows.
    """

    left: str
    right: str
    left_flank: str
    right_flank: str
    unit_bp: int

    @classmethod
    def from_dictionary(
        cls, rd: RepeatDictionary, settings: HybridSettings, layout: ReferenceLayoutSettings
    ) -> Anchors:
        """Anchor on the layout's outer fixed repeats (``left_anchor_id``/``right_anchor_id``)."""
        return cls(
            rd.repeats[layout.left_anchor_id],
            rd.repeats[layout.right_anchor_id],
            rd.flanking_left[-settings.flank_anchor_bp :],
            rd.flanking_right[: settings.flank_anchor_bp],
            rd.repeat_length_bp,
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
    return sum(ord(c) - PHRED_OFFSET for c in qual) / len(qual) if qual else 0.0


def _span_in(
    target: str, anchors: Anchors, settings: HybridSettings
) -> tuple[int, int, int, str] | None:
    """Return (start, end_excl, edits, basis) of motif1..motif9 in target, or None."""
    k = settings.anchor_max_edits
    left = infix_hit(anchors.left, target, k)
    right = None
    if left is not None:
        tail = infix_hit(anchors.right, target[left[1] :], k)
        right = None if tail is None else (tail[0] + left[1], tail[1] + left[1], tail[2])
    if left is not None and right is not None:
        return left[0], right[1], left[2] + right[2], "motif"
    kf = max(settings.flank_anchor_edit_floor, k // settings.flank_anchor_edit_divisor)
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
    lo = settings.min_span_units * anchors.unit_bp
    hi = settings.max_span_units * anchors.unit_bp
    cats = ReadCategories()
    for rec in reads:
        seq = rec.seq.upper()
        best: tuple[tuple[int, int, int, str], str, str, str] | None = None
        for strand, target, qual in (("+", seq, rec.qual), ("-", rc(seq), rec.qual[::-1])):
            hit = _span_in(target, anchors, settings)
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
