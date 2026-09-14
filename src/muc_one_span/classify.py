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

from muc_one_span.classification_summary import (
    _compute_classification_summary as _compute_classification_summary,
)
from muc_one_span.classification_summary import (
    _qual_to_confidence as _qual_to_confidence,
)
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
from muc_one_span.settings import DEFAULT_SETTINGS, ClassificationSettings, ConfidenceSettings

logger = logging.getLogger(__name__)


def classify_repeat(
    sequence: str,
    repeat_dict: RepeatDictionary,
    *,
    settings: ClassificationSettings | None = None,
) -> dict:
    """Classify a single repeat unit against the known dictionary.

    Args:
        sequence: The 60bp (or near-60bp) repeat sequence.
        repeat_dict: The loaded repeat dictionary.

    Returns:
        Classification result dict with type, match, and (for unknowns)
        closest_match, edit_distance, identity_pct, differences.
    """
    effective = settings or DEFAULT_SETTINGS.classification

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

    # Signed net change determines the downstream reading frame.
    net_indel_bases = _compute_net_indel(diffs)
    is_frameshift = has_indels and (net_indel_bases % 3 != 0)

    result: dict = {
        "type": "unknown",
        "match": "closest",
        "closest_match": best_id,
        "edit_distance": best_dist,
        "identity_pct": identity_pct,
        "confidence": identity_pct / 100,
        "differences": diffs,
        "net_indel_bases": net_indel_bases,
    }

    if has_indels:
        result["classification"] = "mutation"
        result["frameshift"] = is_frameshift
    else:
        result["classification"] = (
            "novel_repeat" if best_dist > effective.novel_repeat_edit_distance else "variant"
        )

    return result


def _probe_sizes_generator(
    unit_length: int,
    remaining: int,
    settings: ClassificationSettings,
    *,
    max_indel_probe: int | None = None,
) -> list[int]:
    """Generate probe sizes: canonical first, then small-to-large."""
    probe_limit = settings.max_indel_probe if max_indel_probe is None else max_indel_probe
    sizes = [min(unit_length, remaining)]
    for ps in range(
        max(unit_length - probe_limit, max(1, int(unit_length * settings.minimum_unit_fraction))),
        min(unit_length + probe_limit + 1, remaining + 1),
    ):
        if ps != unit_length:
            sizes.append(ps)
    return sizes


def _classify_backward(
    sequence: str,
    repeat_dict: RepeatDictionary,
    stop_pos: int,
    *,
    settings: ClassificationSettings | None = None,
) -> list[tuple[dict, int, int]]:
    """Classify repeats from 3' end backward, anchored on after-repeats.

    Returns list of (result, start_pos, end_pos) tuples in forward order.
    Stops when reaching stop_pos or when confidence drops.
    """
    effective = settings or DEFAULT_SETTINGS.classification
    unit_length = repeat_dict.repeat_length_bp
    minimum_size = max(1, int(unit_length * effective.minimum_unit_fraction))
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
    while pos - minimum_size > stop_pos:
        remaining_back = pos - stop_pos
        if remaining_back < minimum_size:
            break

        best_result: dict | None = None
        best_size = unit_length
        best_dist: float = float("inf")

        for probe_size in _probe_sizes_generator(unit_length, remaining_back, effective):
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
            best_result = classify_repeat(window, repeat_dict, settings=effective)
            best_size = len(window)
            best_dist = best_result.get("edit_distance", 999)

        if best_result is None or best_dist > effective.max_fit_edit_distance:
            break

        results.append((best_result, pos - best_size, pos))
        pos -= best_size

    results.reverse()
    return results


def _forward_classify(
    sequence: str,
    repeat_dict: RepeatDictionary,
    unit_length: int,
    max_indel_probe: int | None = None,
    strict_segmentation: bool | None = None,
    *,
    settings: ClassificationSettings | None = None,
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
    effective = settings or DEFAULT_SETTINGS.classification
    probe_limit = effective.max_indel_probe if max_indel_probe is None else max_indel_probe
    strict = effective.strict_segmentation if strict_segmentation is None else strict_segmentation
    minimum_size = max(1, int(unit_length * effective.minimum_unit_fraction))
    repeats: list[dict] = []
    mutations: list[dict] = []
    labels: list[str] = []

    pos = 0
    repeat_index = 0
    cumulative_offset = 0

    while pos < len(sequence):
        repeat_index += 1

        remaining = len(sequence) - pos
        if remaining < minimum_size:
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

        for probe_size in _probe_sizes_generator(
            unit_length, remaining, effective, max_indel_probe=probe_limit
        ):
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
                result = classify_repeat(window, repeat_dict, settings=effective)
                if result is not None:
                    dist = 0 if result["match"] == "exact" else result.get("edit_distance", 999)
                    best_dist = dist
                    best_result = result
                    best_window_size = unit_length

            if best_dist > 0:
                for probe_size in range(
                    max(unit_length - probe_limit, minimum_size),
                    min(unit_length + probe_limit + 1, remaining + 1),
                ):
                    if probe_size == unit_length:
                        continue
                    window = sequence[pos : pos + probe_size]
                    result = classify_repeat(window, repeat_dict, settings=effective)
                    if result is None:
                        continue
                    dist = 0 if result["match"] == "exact" else result.get("edit_distance", 999)
                    if dist < best_dist:
                        best_dist = dist
                        best_result = result
                        best_window_size = probe_size
                    if best_dist <= effective.early_stop_edit_distance:
                        break

        if best_result is None:
            raise RuntimeError(
                f"Classification failed: no match found at position {pos} "
                f"(remaining: {remaining} bp, repeat index: {repeat_index})"
            )
        if strict and (
            best_dist > effective.max_fit_edit_distance
            or any(b not in "ACGT" for b in sequence[pos : pos + best_window_size])
        ):
            break
        result = best_result
        advance = best_window_size

        # Track cumulative offset if this window is non-standard size
        if best_window_size != unit_length:
            net_indel = best_window_size - unit_length
            cumulative_offset += net_indel

        result["index"] = repeat_index
        result["start"] = pos
        result["end"] = pos + advance
        result["fit_status"] = (
            "unresolved"
            if best_dist > effective.max_fit_edit_distance
            or any(b not in "ACGT" for b in sequence[pos : pos + advance])
            else "accepted"
        )

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
                        "frameshift": (
                            best_window_size - len(repeat_dict.repeats[result["parent_repeat"]])
                        )
                        % 3
                        != 0,
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
    # Kept as a compatibility import. An unknown gap cannot safely be treated as
    # one mutated repeat or assigned a biological repeat count.
    return repeats, mutations, labels


def classify_sequence(
    sequence: str,
    repeat_dict: RepeatDictionary,
    *,
    strict_segmentation: bool | None = None,
    settings: ClassificationSettings | None = None,
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
    effective = settings or DEFAULT_SETTINGS.classification
    strict = effective.strict_segmentation if strict_segmentation is None else strict_segmentation
    unit_length = repeat_dict.repeat_length_bp

    repeats, mutations, labels, pos, cumulative_offset = _forward_classify(
        sequence, repeat_dict, unit_length, strict_segmentation=strict, settings=effective
    )

    # A variable-length template immediately before an unresolved gap can borrow
    # the gap's first bases (e.g. X followed by a long A insertion looks like dupA).
    # Keep it as an unresolved candidate, not a localized mutation assertion.
    unresolved_candidates = []
    if strict and pos < len(sequence) and repeats and repeats[-1].get("mutation_name"):
        last = repeats.pop()
        labels.pop()
        unresolved_candidates = [m for m in mutations if m["repeat_index"] == last["index"]]
        mutations = [m for m in mutations if m["repeat_index"] != last["index"]]
        cumulative_offset -= last["end"] - last["start"] - unit_length
        pos = last["start"]
    result = _compute_classification_summary(repeats, mutations, labels, cumulative_offset)
    suffix: list[dict] = []
    suffix_start = len(sequence)
    if pos < len(sequence):
        backward = _classify_backward(sequence, repeat_dict, pos, settings=effective)
        for repeat, start, end in reversed(backward):
            if repeat["match"] != "exact" or end != suffix_start:
                break
            suffix.append({**repeat, "start": start, "end": end, "index": None})
            suffix_start = start
        suffix.reverse()
    unresolved = (
        [{"start": pos, "end": suffix_start, "reason": "unresolved_sequence"}]
        if pos < suffix_start
        else []
    )
    unresolved_fit = [r for r in repeats if r.get("fit_status") == "unresolved"]
    uncertain_index = min((r["index"] for r in unresolved_fit), default=len(repeats) + 1)
    for repeat in repeats:
        repeat["localization_status"] = (
            "ambiguous" if repeat["index"] >= uncertain_index else "resolved"
        )
    for mutation in mutations:
        mutation["localization_status"] = (
            "ambiguous" if mutation["repeat_index"] >= uncertain_index else "resolved"
        )
    result.update(
        {
            "sequence_length": len(sequence),
            "classified_bases": pos,
            "classification_coverage": pos / len(sequence) if sequence else 0.0,
            "unclassified_regions": unresolved,
            "unresolved_candidates": unresolved_candidates,
            "unresolved_regions": [
                *unresolved,
                *[
                    {"start": r["start"], "end": r["end"], "reason": "uncertain_repeat_fit"}
                    for r in unresolved_fit
                ],
            ],
            "segmentation_policy": "strict_experimental" if strict else "candidate_windows",
            "recovered_suffix": suffix,
            "ambiguous_bases": sum(b not in "ACGT" for b in sequence),
            "reconstruction_status": (
                "complete_segmentation"
                if pos == len(sequence) and sequence and not unresolved_fit
                else "ambiguous_reconstruction"
                if sequence
                else "insufficient_evidence"
            ),
            "confidence_semantics": "heuristic_dictionary_fit_not_probability",
        }
    )
    return result


def validate_mutations_against_vcf(
    classification_result: dict,
    vcf_variants: list[dict] | None = None,
    flank_length: int = 500,
    unit_length: int = 60,
    boundary_repeats: int | None = None,
    boundary_penalty: float | None = None,
    *,
    sequence: str | None = None,
    repeat_dict: RepeatDictionary | None = None,
    consensus_context: dict | None = None,
    settings: ConfidenceSettings | None = None,
) -> dict:
    """Annotate exact VCF sequence concordance and heuristic evidence weights.

    Exact support requires actual reference/consensus replay context, the observed
    VNTR sequence and dictionary. Legacy position/QUAL-only records cannot prove
    variant identity. This is concordance with the source VCF, not independent
    validation. ``flank_length``/``unit_length`` remain compatibility arguments.
    """
    from muc_one_span.variant_support import mutation_concordance

    result = classification_result.copy()
    result["mutations_detected"] = [m.copy() for m in result.get("mutations_detected", [])]
    result["repeats"] = [r.copy() for r in result.get("repeats", [])]

    if vcf_variants is None:
        return result

    effective = settings or DEFAULT_SETTINGS.confidence
    boundary_count = effective.boundary_repeats if boundary_repeats is None else boundary_repeats
    boundary_weight = effective.boundary_penalty if boundary_penalty is None else boundary_penalty
    total_repeats = len(result["repeats"])
    support = mutation_concordance(result, vcf_variants, sequence, repeat_dict, consensus_context)

    for mutation in result["mutations_detected"]:
        repeat_idx = mutation["repeat_index"]
        status, supporting = support[repeat_idx]
        mutation["vcf_support"] = bool(supporting)
        mutation["vcf_support_status"] = status
        mutation["vcf_qual"] = max((v.get("qual") or 0.0 for v in supporting), default=0.0)

        # Check if mutation is near the allele boundary.
        # Only apply to alleles long enough for boundary to be meaningful
        # (at least 2x boundary_repeats).
        is_boundary = (
            total_repeats > 2 * boundary_count and repeat_idx > total_repeats - boundary_count
        )
        mutation["boundary"] = is_boundary

        # Adjust confidence in the corresponding repeat
        if repeat_idx - 1 < len(result["repeats"]):
            repeat_result = result["repeats"][repeat_idx - 1]
            base_confidence = repeat_result.get("confidence", 1.0)
            if supporting:
                vcf_score = _qual_to_confidence(mutation["vcf_qual"], settings=effective)
            elif status == "absent":
                vcf_score = effective.absent_weight
            else:
                vcf_score = 1.0
            confidence = base_confidence * vcf_score
            # Apply boundary penalty for mutations near allele ends
            if is_boundary:
                confidence *= boundary_weight
            repeat_result["confidence"] = round(confidence, 4)

    # Recompute allele_confidence
    confidences = [r.get("confidence", 1.0) for r in result["repeats"]]
    result["allele_confidence"] = (
        round(sum(confidences) / len(confidences), 4) if confidences else 0.0
    )

    return result
