"""Diagnostic terminal-anchor spans; these do not replace allele-length calling.

Input record indices identify observations, even when read names collide. They
are not proof of independent PCR molecules. Supply source records once, excluding
secondary/supplementary alignment copies, and retain any source identity sidecar.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from muc_one_span.repeat_alignment import edit_distance

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class SpanEvidence:
    """Observed span interval in oriented-read bases, not calibrated uncertainty.

    The interval covers equally optimal anchor boundaries at one local locus.
    It does not correct interior sequencing indels or establish repeat count.
    Distinct anchor loci and orientations are rejected as ambiguous.
    """

    record_index: int
    read_name: str
    status: str
    orientation: str | None = None
    span_min: int | None = None
    span_max: int | None = None
    anchor_edits: int | None = None


def _anchor_hits(sequence: str, anchor: str, max_edits: int) -> list[tuple[int, int, int]]:
    """Return all minimum-edit substring hits, including tied boundary choices."""
    exact: list[tuple[int, int, int]] = []
    start = sequence.find(anchor)
    while start >= 0:
        exact.append((start, start + len(anchor), 0))
        start = sequence.find(anchor, start + 1)
    if exact or max_edits == 0:
        return exact

    # Myers recurrence with a free text prefix (no low-bit insertion on PH).
    masks: dict[str, int] = {}
    for index, symbol in enumerate(anchor):
        masks[symbol] = masks.get(symbol, 0) | (1 << index)
    positive, negative, distance = ~0, 0, len(anchor)
    final_bit = 1 << (len(anchor) - 1)
    best = max_edits
    hits: list[tuple[int, int, int]] = []
    for end, symbol in enumerate(sequence, 1):
        matches = masks.get(symbol, 0)
        vertical = matches | negative
        horizontal = (((matches & positive) + positive) ^ positive) | matches
        ph = negative | ~(horizontal | positive)
        mh = positive & horizontal
        distance += bool(ph & final_bit) - bool(mh & final_bit)
        ph <<= 1
        mh <<= 1
        positive = mh | ~(vertical | ph)
        negative = ph & vertical
        if distance > best:
            continue
        if distance < best:
            best, hits = distance, []
        for width in range(max(1, len(anchor) - best), len(anchor) + best + 1):
            start = end - width
            if start >= 0 and edit_distance(anchor, sequence[start:end]) == best:
                hits.append((start, end, best))
    return hits


def spanning_evidence(
    records: Iterable[tuple[str, str]],
    left_anchor: str,
    right_anchor: str,
    *,
    max_edits: int = 1,
) -> list[SpanEvidence]:
    """Scan each input record in both orientations for diagnostic spanning evidence.

    Anchors must be nonempty with a per-anchor edit allowance smaller than either
    anchor. Matching uses literal Levenshtein distance, including substitutions
    and indels. Only globally best matches for each anchor are retained; this
    conservative rule can reject reads with a degraded true anchor and a better
    ectopic match. Terminal-anchor mutations and unrecognized partial reads can
    cause ascertainment bias. No truth sequences or allele labels are consulted.
    """
    if (
        not left_anchor
        or not right_anchor
        or not 0 <= max_edits < min(len(left_anchor), len(right_anchor))
    ):
        raise ValueError("Require nonempty anchors and a smaller nonnegative edit allowance")
    results: list[SpanEvidence] = []
    for record_index, (name, sequence) in enumerate(records):
        pairs: list[tuple[str, int, int, int]] = []
        partial = False
        for orientation, oriented in (
            ("forward", sequence),
            ("reverse", sequence.translate(_COMPLEMENT)[::-1]),
        ):
            left = _anchor_hits(oriented, left_anchor, max_edits)
            right = _anchor_hits(oriented, right_anchor, max_edits)
            partial |= bool(left or right)
            pairs.extend(
                (orientation, first[0], last[1], first[2] + last[2])
                for first in left
                for last in right
                if first[1] <= last[0]
            )
        if not pairs:
            status = "partial_anchor_pair" if partial else "no_anchor_pair"
            results.append(SpanEvidence(record_index, name, status))
            continue
        starts, ends = [p[1] for p in pairs], [p[2] for p in pairs]
        if (
            len({p[0] for p in pairs}) > 1
            or max(starts) - min(starts) > 2 * max_edits
            or max(ends) - min(ends) > 2 * max_edits
        ):
            results.append(SpanEvidence(record_index, name, "ambiguous_anchor_pair"))
            continue
        spans = [p[2] - p[1] for p in pairs]
        results.append(
            SpanEvidence(
                record_index, name, "spanning", pairs[0][0], min(spans), max(spans), pairs[0][3]
            )
        )
    return results


def assign_span(
    evidence: SpanEvidence, candidate_spans: dict[str, int], max_distance: int
) -> tuple[str, ...]:
    """Return every length-compatible candidate; never break a haplotype tie.

    Candidate spans are base lengths. All candidates intersecting the observed
    boundary interval expanded by ``max_distance`` are returned, not just the
    nearest. Equal length cannot establish sequence identity or assignment purity.
    This is diagnostic compatibility, not a validated haplotype assignment model.
    """
    if max_distance < 0:
        raise ValueError("max_distance must be nonnegative")
    if evidence.status != "spanning" or evidence.span_min is None or evidence.span_max is None:
        return ()
    return tuple(
        name
        for name, span in candidate_spans.items()
        if evidence.span_min - max_distance <= span <= evidence.span_max + max_distance
    )
