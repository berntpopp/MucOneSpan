"""S5/S6: ladder-flanked allele references and edit-distance competition for every read.

``assign_flank_bp``/``assign_margin``/``assign_max_error_rate``/``min_fragment_bp`` are
validated ``HybridSettings`` fields (``muc_one_span.settings``); nothing here hardcodes a
threshold and no parameter has a default, so callers always pass the values of the
settings they run with.
"""

from __future__ import annotations

from dataclasses import dataclass

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import edit_distance_infix, project, rc
from muc_one_span.hybrid.spans import ReadRecord
from muc_one_span.settings import HybridSettings

OFF_TARGET = "off_target"
UNDECIDED = "undecided"


def hybrid_references(
    drafts: dict[str, str],
    rd: RepeatDictionary,
    flank_bp: int,
) -> dict[str, str]:
    """Wrap each motif1..motif9 draft in the ladder's hg38 flanks.

    ``flank_bp`` has no default: callers pass ``HybridSettings.assign_flank_bp``.
    """
    left, right = rd.flanking_left[-flank_bp:], rd.flanking_right[:flank_bp]
    return {name: left + seq + right for name, seq in drafts.items()}


@dataclass(frozen=True)
class Assignment:
    """Outcome for one read: an allele name, ``undecided``, or ``off_target``."""

    allele: str | None
    margin: int
    distances: dict[str, int]
    oriented: str


def assign_read(seq: str, refs: dict[str, str], margin: int, max_error_rate: float) -> Assignment:
    """Assign seq to the reference with the smallest infix edit distance, if clearly best.

    Returns ``allele=None`` when the margin is too small, and ``allele="off_target"`` when
    even the best reference needs more than ``max_error_rate * len(seq)`` edits. A single
    reference (homozygous sample) is assigned subject to the same off-target guard.

    ``margin`` and ``max_error_rate`` have no default: callers source them from
    ``HybridSettings.assign_margin``/``assign_max_error_rate`` (never a hidden literal).
    """
    best: Assignment | None = None
    for oriented in (seq, rc(seq)):
        dist = {name: edit_distance_infix(oriented, ref) for name, ref in refs.items()}
        ranked = sorted(dist.values())
        gap = ranked[1] - ranked[0] if len(ranked) > 1 else len(seq)
        winner = min(dist, key=lambda k: dist[k])
        cand = Assignment(winner if gap >= margin else None, gap, dist, oriented)
        if best is None or ranked[0] < min(best.distances.values()):
            best = cand
    if best is None:
        raise ValueError("no reference to assign against")
    if min(best.distances.values()) > max_error_rate * len(seq):
        return Assignment(OFF_TARGET, best.margin, best.distances, best.oriented)
    return best


def assign_reads(
    reads: list[ReadRecord], refs: dict[str, str], settings: HybridSettings
) -> dict[str, list[str]]:
    """Map allele -> oriented read sequences, plus ``undecided`` and ``off_target``.

    Reads shorter than ``settings.min_fragment_bp`` are not considered.
    """
    out: dict[str, list[str]] = {name: [] for name in refs}
    out[UNDECIDED] = []
    out[OFF_TARGET] = []
    for rec in reads:
        if len(rec.seq) < settings.min_fragment_bp:
            continue
        got = assign_read(
            rec.seq.upper(), refs, settings.assign_margin, settings.assign_max_error_rate
        )
        out[got.allele or UNDECIDED].append(got.oriented)
    return out


def trim_to_draft(oriented: str, ref: str, flank_bp: int, draft_len: int) -> str:
    """Cut an assigned read to the part aligned inside the draft (drop ladder flanks).

    ``flank_bp`` has no default: callers pass the same value used to build ``ref`` via
    ``hybrid_references``, sourced from ``HybridSettings.assign_flank_bp``.
    """
    proj = project(oriented, ref, partial=True)
    lo, hi = max(proj.first, flank_bp), min(proj.last, flank_bp + draft_len)
    if hi - lo <= 0:
        return ""
    return oriented[proj.t2q[lo] : proj.t2q[hi]]
