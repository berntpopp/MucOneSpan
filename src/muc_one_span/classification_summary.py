"""Legacy dictionary-fit summaries and heuristic VCF weights (not probabilities)."""

from __future__ import annotations

from muc_one_span.settings import DEFAULT_SETTINGS, ConfidenceSettings


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


def _qual_to_confidence(qual: float, *, settings: ConfidenceSettings | None = None) -> float:
    """Map VCF QUAL score to a confidence weight in [0, 1].

    Uses configured low/high breakpoints and linear interpolation between their
    configured weights. Values below the low breakpoint use ``weight_below``.

    Args:
        qual: VCF QUAL score.
        settings: Immutable confidence settings; central defaults when omitted.

    Returns:
        Confidence weight between 0.3 and 1.0.
    """
    effective = settings or DEFAULT_SETTINGS.confidence
    if qual >= effective.qual_high:
        return float(effective.weight_high)
    if qual >= effective.qual_low:
        fraction = (qual - effective.qual_low) / (effective.qual_high - effective.qual_low)
        return effective.weight_low + fraction * (effective.weight_high - effective.weight_low)
    return effective.weight_below
