# src/muc_one_span/classify.py
"""Repeat unit classification and mutation detection.

The classification algorithm handles frameshifted sequences by tracking
the cumulative indel offset.  When a repeat contains an insertion or
deletion, subsequent window boundaries are shifted by the net indel
length so that downstream repeats are correctly framed.

For example, a dupC (1bp insertion) at repeat 25 shifts all windows
after repeat 25 by +1bp.  Without correction, every downstream window
would straddle two repeat boundaries and fail to match any known type.
With correction, the windows realign to the true repeat boundaries and
classify correctly.
"""

from __future__ import annotations

import logging

from muc_one_span.classify_types import (
    MutationDetected as MutationDetected,
)
from muc_one_span.classify_types import (
    RepeatClassification as RepeatClassification,
)
from muc_one_span.classify_types import (
    RepeatDifference as RepeatDifference,
)
from muc_one_span.classify_types import (
    SequenceClassification as SequenceClassification,
)
from muc_one_span.config import RepeatDictionary
from muc_one_span.repeat_alignment import (
    _compute_net_indel as _compute_net_indel,
)
from muc_one_span.repeat_alignment import (
    characterize_differences as characterize_differences,
)
from muc_one_span.repeat_alignment import (
    edit_distance as edit_distance,
)

logger = logging.getLogger(__name__)


def classify_repeat(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> dict:
    """Classify a single repeat unit against the known dictionary.

    Args:
        sequence: The 60bp (or near-60bp) repeat sequence.
        repeat_dict: The loaded repeat dictionary.

    Returns:
        Classification result dict with type, match, and (for unknowns)
        closest_match, edit_distance, identity_pct, differences.
    """
    # O(1) exact-match lookup via cached reverse map (sequence -> ID)
    if sequence in repeat_dict.seq_to_id:
        return {"type": repeat_dict.seq_to_id[sequence], "match": "exact", "confidence": 1.0}

    # Check mutation templates (variable-length exact matches)
    if repeat_dict.mutated_sequences and sequence in repeat_dict.mutated_sequences:
        parent_repeat, mut_name = repeat_dict.mutated_sequences[sequence]
        return {
            "type": f"{parent_repeat}:{mut_name}",
            "match": "exact",
            "confidence": 1.0,
            "mutation_name": mut_name,
            "parent_repeat": parent_repeat,
        }

    # No exact match -- find closest by edit distance
    best_id = ""
    best_dist = float("inf")

    for repeat_id, ref_seq in repeat_dict.repeats.items():
        dist = edit_distance(ref_seq, sequence)
        if dist < best_dist:
            best_dist = dist
            best_id = repeat_id

    # Characterize the specific differences
    ref_seq = repeat_dict.repeats[best_id]
    diffs = characterize_differences(ref_seq, sequence)

    # Calculate identity percentage based on alignment length
    max_len = max(len(ref_seq), len(sequence))
    identity_pct = round((1 - best_dist / max_len) * 100, 1) if max_len > 0 else 0.0

    # Determine if differences contain indels
    has_indels = any(d["type"] in ("insertion", "deletion") for d in diffs)

    # Calculate total indel length for frameshift check
    indel_bases = sum(
        len(d["alt"]) if d["type"] == "insertion" else len(d["ref"])
        for d in diffs
        if d["type"] in ("insertion", "deletion")
    )
    is_frameshift = has_indels and (indel_bases % 3 != 0)

    result: dict = {
        "type": "unknown",
        "match": "closest",
        "closest_match": best_id,
        "edit_distance": best_dist,
        "identity_pct": identity_pct,
        "confidence": identity_pct / 100,
        "differences": diffs,
    }

    if has_indels:
        result["classification"] = "mutation"
        result["frameshift"] = is_frameshift
    else:
        result["classification"] = "novel_repeat" if best_dist > 2 else "variant"

    return result


def _probe_sizes_generator(unit_length: int, max_indel_probe: int, remaining: int) -> list[int]:
    """Generate probe sizes: canonical first, then small-to-large."""
    sizes = [min(unit_length, remaining)]
    for ps in range(
        max(unit_length - max_indel_probe, unit_length // 2),
        min(unit_length + max_indel_probe + 1, remaining + 1),
    ):
        if ps != unit_length:
            sizes.append(ps)
    return sizes


def _classify_backward(
    sequence: str,
    repeat_dict: RepeatDictionary,
    stop_pos: int,
) -> list[tuple[dict, int, int]]:
    """Classify repeats from 3' end backward, anchored on after-repeats.

    Returns list of (result, start_pos, end_pos) tuples in forward order.
    Stops when reaching stop_pos or when confidence drops.
    """
    unit_length = repeat_dict.repeat_length_bp
    max_indel_probe = 30
    after_ids = list(reversed(repeat_dict.after_repeat_ids))

    results: list[tuple[dict, int, int]] = []
    pos = len(sequence)

    # First try to match after-repeats from the end
    for expected_id in after_ids:
        if pos - unit_length < stop_pos:
            break
        window = sequence[pos - unit_length : pos]
        if window in repeat_dict.seq_to_id and repeat_dict.seq_to_id[window] == expected_id:
            result = {"type": expected_id, "match": "exact", "confidence": 1.0}
            results.append((result, pos - unit_length, pos))
            pos -= unit_length
        else:
            break

    # Continue backward through canonical region
    while pos - unit_length // 2 > stop_pos:
        remaining_back = pos - stop_pos
        if remaining_back < unit_length // 2:
            break

        best_result: dict | None = None
        best_size = unit_length
        best_dist: float = float("inf")

        for probe_size in _probe_sizes_generator(unit_length, max_indel_probe, remaining_back):
            start = pos - probe_size
            if start < stop_pos:
                continue
            window = sequence[start:pos]
            if window in repeat_dict.seq_to_id:
                best_result = {
                    "type": repeat_dict.seq_to_id[window],
                    "match": "exact",
                    "confidence": 1.0,
                }
                best_size = probe_size
                best_dist = 0
                break
            if repeat_dict.mutated_sequences and window in repeat_dict.mutated_sequences:
                parent, mname = repeat_dict.mutated_sequences[window]
                best_result = {
                    "type": f"{parent}:{mname}",
                    "match": "exact",
                    "confidence": 1.0,
                    "mutation_name": mname,
                    "parent_repeat": parent,
                }
                best_size = probe_size
                best_dist = 0
                break

        if best_dist > 0:
            # Edit distance fallback
            window = sequence[max(stop_pos, pos - unit_length) : pos]
            best_result = classify_repeat(window, repeat_dict)
            best_size = len(window)
            best_dist = best_result.get("edit_distance", 999)

        if best_result is None or best_dist > 3:
            break

        results.append((best_result, pos - best_size, pos))
        pos -= best_size

    results.reverse()
    return results


def _forward_classify(
    sequence: str,
    repeat_dict: RepeatDictionary,
    unit_length: int,
    max_indel_probe: int,
) -> tuple[list[dict], list[dict], list[str], int, int]:
    """Classify repeats in a forward pass from 5' to 3'.

    Args:
        sequence: Full consensus sequence.
        repeat_dict: The loaded repeat dictionary.
        unit_length: Expected repeat unit length in bp.
        max_indel_probe: Maximum indel size to probe on either side.

    Returns:
        Tuple of (repeats, mutations, labels, pos, cumulative_offset) where
        *pos* is the position where the forward pass stopped and
        *cumulative_offset* is the total net indel accumulated.
    """
    repeats: list[dict] = []
    mutations: list[dict] = []
    labels: list[str] = []

    pos = 0
    repeat_index = 0
    cumulative_offset = 0

    while pos < len(sequence):
        repeat_index += 1

        remaining = len(sequence) - pos
        if remaining < unit_length // 2:
            break

        best_result: dict | None = None
        best_dist = float("inf")
        best_window_size = unit_length

        # --- Phase 1: Check ALL probe sizes for exact match ---
        # First check canonical size for standard repeats (common case).
        # Then check ALL sizes for mutation templates (which are non-60bp).
        # Mutation templates take priority over canonical-size standard matches
        # because they explain the actual biological repeat length.
        exact_found = False
        canonical_result: dict | None = None

        for probe_size in _probe_sizes_generator(unit_length, max_indel_probe, remaining):
            window = sequence[pos : pos + probe_size]
            # Check mutation templates first (variable-length exact matches)
            if repeat_dict.mutated_sequences and window in repeat_dict.mutated_sequences:
                parent, mname = repeat_dict.mutated_sequences[window]
                best_result = {
                    "type": f"{parent}:{mname}",
                    "match": "exact",
                    "confidence": 1.0,
                    "mutation_name": mname,
                    "parent_repeat": parent,
                }
                best_window_size = probe_size
                best_dist = 0
                exact_found = True
                break
            # Check standard repeats
            if window in repeat_dict.seq_to_id and canonical_result is None:
                canonical_result = {
                    "type": repeat_dict.seq_to_id[window],
                    "match": "exact",
                    "confidence": 1.0,
                }

        # Use canonical match if no mutation template found
        if not exact_found and canonical_result is not None:
            best_result = canonical_result
            best_window_size = unit_length
            best_dist = 0
            exact_found = True

        # --- Phase 2: Edit distance fallback (only if no exact match) ---
        if not exact_found:
            # Try canonical size first
            if remaining >= unit_length:
                window = sequence[pos : pos + unit_length]
                result = classify_repeat(window, repeat_dict)
                if result is not None:
                    dist = 0 if result["match"] == "exact" else result.get("edit_distance", 999)
                    best_dist = dist
                    best_result = result
                    best_window_size = unit_length

            if best_dist > 0:
                for probe_size in range(
                    max(unit_length - max_indel_probe, unit_length // 2),
                    min(unit_length + max_indel_probe + 1, remaining + 1),
                ):
                    if probe_size == unit_length:
                        continue
                    window = sequence[pos : pos + probe_size]
                    result = classify_repeat(window, repeat_dict)
                    if result is None:
                        continue
                    dist = 0 if result["match"] == "exact" else result.get("edit_distance", 999)
                    if dist < best_dist:
                        best_dist = dist
                        best_result = result
                        best_window_size = probe_size
                    if best_dist <= 1:
                        break

        if best_result is None:
            raise RuntimeError(
                f"Classification failed: no match found at position {pos} "
                f"(remaining: {remaining} bp, repeat index: {repeat_index})"
            )
        result = best_result
        advance = best_window_size

        # Track cumulative offset if this window is non-standard size
        if best_window_size != unit_length:
            net_indel = best_window_size - unit_length
            cumulative_offset += net_indel

        result["index"] = repeat_index

        if result["match"] == "exact":
            labels.append(result["type"])
            # Track template-matched mutations in mutations_detected
            if result.get("mutation_name"):
                mutations.append(
                    {
                        "repeat_index": repeat_index,
                        "closest_type": result.get("parent_repeat", result["type"]),
                        "mutation_name": result["mutation_name"],
                        "template_match": True,
                        "frameshift": True,  # known mutations are frameshifts
                    }
                )
        elif result.get("classification") == "mutation":
            # Use MucOneUp nomenclature: "Xm" = repeat X with mutation
            labels.append(f"{result['closest_match']}m")
            mutations.append(
                {
                    "repeat_index": repeat_index,
                    "closest_type": result["closest_match"],
                    "differences": result["differences"],
                    "frameshift": result.get("frameshift", False),
                }
            )
        else:
            # Variant of known type (substitutions only, no indel)
            labels.append(f"?{result.get('closest_match', '?')}")

        repeats.append(result)
        pos += max(advance, 1)  # always advance at least 1 to avoid infinite loop

    return repeats, mutations, labels, pos, cumulative_offset


def _apply_bidirectional_fallback(
    sequence: str,
    repeat_dict: RepeatDictionary,
    repeats: list[dict],
    mutations: list[dict],
    labels: list[str],
    forward_pos: int,
) -> tuple[list[dict], list[dict], list[str]]:
    """Apply bidirectional fallback when the forward pass left unconsumed sequence.

    If the forward pass stopped with significant unconsumed sequence,
    classify from the 3' end backward and bridge the gap.

    Args:
        sequence: Full consensus sequence.
        repeat_dict: The loaded repeat dictionary.
        repeats: Repeat classifications accumulated by the forward pass (mutated in place).
        mutations: Mutations accumulated by the forward pass (mutated in place).
        labels: Labels accumulated by the forward pass (mutated in place).
        forward_pos: Position where the forward pass stopped.

    Returns:
        Tuple of (repeats, mutations, labels) with fallback results appended.
    """
    unit_length = repeat_dict.repeat_length_bp
    pos = forward_pos
    repeat_index = len(repeats)

    # --- Bidirectional fallback ---
    # If the forward pass stopped with significant unconsumed sequence,
    # try classifying from the 3' end backward.
    if pos < len(sequence) - unit_length // 2:
        backward = _classify_backward(sequence, repeat_dict, pos)
        if backward:
            # Gap between forward and backward = mutated region
            gap_start = pos
            gap_end = backward[0][1]
            if gap_end > gap_start:
                gap_seq = sequence[gap_start:gap_end]
                gap_result = classify_repeat(gap_seq, repeat_dict)
                repeat_index += 1
                gap_result["index"] = repeat_index
                if gap_result.get("classification") == "mutation" or gap_result["match"] != "exact":
                    labels.append(f"{gap_result.get('closest_match', '?')}m")
                    mutations.append(
                        {
                            "repeat_index": repeat_index,
                            "closest_type": gap_result.get("closest_match", "?"),
                            "differences": gap_result.get("differences", []),
                            "frameshift": gap_result.get("frameshift", False),
                        }
                    )
                else:
                    labels.append(gap_result["type"])
                repeats.append(gap_result)

            # Append backward results
            for bwd_result, _bwd_start, _bwd_end in backward:
                repeat_index += 1
                bwd_result["index"] = repeat_index
                if bwd_result["match"] == "exact":
                    labels.append(bwd_result["type"])
                else:
                    labels.append(f"?{bwd_result.get('closest_match', '?')}")
                repeats.append(bwd_result)

    return repeats, mutations, labels


def _compute_classification_summary(
    repeats: list[dict],
    mutations: list[dict],
    labels: list[str],
    cumulative_offset: int,
) -> dict:
    """Compute summary statistics and build the final classification result dict.

    Args:
        repeats: Per-repeat classification results.
        mutations: Mutations detected during classification.
        labels: Label string for each repeat.
        cumulative_offset: Net cumulative indel offset accumulated during classification.

    Returns:
        Final classification result dict.
    """
    confidences = [r.get("confidence", 1.0) for r in repeats]
    exact_count = sum(1 for r in repeats if r.get("match") == "exact")
    allele_confidence = sum(confidences) / len(confidences) if confidences else 0.0
    exact_match_pct = (exact_count / len(repeats) * 100) if repeats else 0.0

    return {
        "structure": " ".join(labels),
        "repeats": repeats,
        "mutations_detected": mutations,
        "cumulative_offset": cumulative_offset,
        "allele_confidence": round(allele_confidence, 4),
        "exact_match_pct": round(exact_match_pct, 1),
    }


def classify_sequence(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> dict:
    """Classify all repeat units in a consensus sequence.

    Uses offset-aware windowing: when a repeat contains an indel, the
    cumulative offset is tracked and subsequent window boundaries are
    shifted accordingly.  This corrects for frameshift propagation --
    a 1bp insertion at repeat 25 would otherwise misalign all downstream
    windows.

    Algorithm:
        1. Start at position 0 with offset = 0
        2. Extract window of ``unit_length + offset`` bases (the mutated
           repeat is longer/shorter than 60bp)
        3. Classify the window
        4. If classification finds indels, compute the net offset and
           accumulate it for subsequent windows
        5. Advance position by ``unit_length + net_indel`` (actual length
           of the repeat in the sequence)
        6. Reset offset to 0 for the next window (each downstream repeat
           is expected to be 60bp again, just starting from the shifted
           position)

    Args:
        sequence: Full consensus sequence (flanking regions should be trimmed).
        repeat_dict: The loaded repeat dictionary.

    Returns:
        Dict with structure string, per-repeat details, and mutation report.
    """
    logger.info("Classifying sequence of %d bp", len(sequence))
    unit_length = repeat_dict.repeat_length_bp
    # Maximum indel size to probe.  Covers all known MUC1 mutations
    # (largest known: 25bp insertion, 14bp deletion).
    max_indel_probe = 30

    repeats, mutations, labels, pos, cumulative_offset = _forward_classify(
        sequence, repeat_dict, unit_length, max_indel_probe
    )

    repeats, mutations, labels = _apply_bidirectional_fallback(
        sequence, repeat_dict, repeats, mutations, labels, pos
    )

    return _compute_classification_summary(repeats, mutations, labels, cumulative_offset)


def _qual_to_confidence(qual: float) -> float:
    """Map VCF QUAL score to a confidence weight in [0, 1].

    Uses linear interpolation between QUAL=5 (0.5) and QUAL=20 (1.0).
    Below QUAL=5, returns 0.3 as a floor. Above 20, returns 1.0.

    This replaces the previous binary threshold (QUAL>=20 → 1.0, else 0.7)
    to give a continuous signal that better reflects Clair3's confidence.

    Args:
        qual: VCF QUAL score.

    Returns:
        Confidence weight between 0.3 and 1.0.
    """
    if qual >= 20.0:
        return 1.0
    if qual >= 5.0:
        return 0.5 + 0.5 * (qual - 5.0) / 15.0
    return 0.3


def validate_mutations_against_vcf(
    classification_result: dict,
    vcf_variants: list[dict] | None = None,
    flank_length: int = 500,
    unit_length: int = 60,
    boundary_repeats: int = 3,
    boundary_penalty: float = 0.5,
) -> dict:
    """Cross-validate detected mutations against VCF variant positions.

    For each mutation, check if a VCF variant overlaps the repeat's
    genomic position.  Adds ``vcf_support`` flag and adjusts confidence
    using the continuous :func:`_qual_to_confidence` function.

    Mutations near the ends of an allele (within *boundary_repeats* of
    the last repeat) receive an additional confidence penalty because
    Clair3 produces systematic artifacts at contig boundaries where
    read alignment quality degrades.

    Args:
        classification_result: Output from :func:`classify_sequence`.
        vcf_variants: List of dicts with ``pos`` (int) and ``qual`` (float).
            If None, VCF validation is skipped (standalone classify mode).
        flank_length: Flanking bp on each side of the contig.
        unit_length: Expected repeat unit length (60bp).
        boundary_repeats: Number of repeats at allele ends subject to
            boundary penalty (default 3).
        boundary_penalty: Confidence multiplier for boundary mutations
            (default 0.5).

    Returns:
        Updated classification result with VCF validation annotations.
    """
    result = classification_result.copy()
    result["mutations_detected"] = [m.copy() for m in result.get("mutations_detected", [])]
    result["repeats"] = [r.copy() for r in result.get("repeats", [])]

    if vcf_variants is None:
        return result

    total_repeats = len(result["repeats"])

    for mutation in result["mutations_detected"]:
        repeat_idx = mutation["repeat_index"]
        # Map repeat index to contig coordinates
        repeat_start = flank_length + (repeat_idx - 1) * unit_length
        repeat_end = repeat_start + unit_length + 30  # allow for indels

        # Check if any VCF variant overlaps this repeat
        supporting = [v for v in vcf_variants if repeat_start <= v["pos"] <= repeat_end]
        mutation["vcf_support"] = len(supporting) > 0
        mutation["vcf_qual"] = max((v["qual"] for v in supporting), default=0.0)

        # Check if mutation is near the allele boundary.
        # Only apply to alleles long enough for boundary to be meaningful
        # (at least 2x boundary_repeats).
        is_boundary = (
            total_repeats > 2 * boundary_repeats and repeat_idx > total_repeats - boundary_repeats
        )
        mutation["boundary"] = is_boundary

        # Adjust confidence in the corresponding repeat
        if repeat_idx - 1 < len(result["repeats"]):
            repeat_result = result["repeats"][repeat_idx - 1]
            base_confidence = repeat_result.get("confidence", 1.0)
            vcf_score = _qual_to_confidence(mutation["vcf_qual"]) if supporting else 0.3
            confidence = base_confidence * vcf_score
            # Apply boundary penalty for mutations near allele ends
            if is_boundary:
                confidence *= boundary_penalty
            repeat_result["confidence"] = round(confidence, 4)

    # Recompute allele_confidence
    confidences = [r.get("confidence", 1.0) for r in result["repeats"]]
    result["allele_confidence"] = (
        round(sum(confidences) / len(confidences), 4) if confidences else 0.0
    )

    return result
