"""Allele-selection and depth gates computed from ladder alignment evidence.

Primary alignment records per contig approximate molecules (one primary record
per read). A cluster whose primary records form a second mode at least
``min_gap`` units from its peak may hide another allele or a fragment class, so
its selection is unresolved. These gates only prevent a reassuring negative;
they never create a variant call.
"""

from __future__ import annotations

from typing import Any

from muc_one_span.settings import DEFAULT_SETTINGS, AlleleSelectionSettings


def _contig_units(name: str) -> int | None:
    """Repeat units of a ladder ``contig_N`` name; None for any other name."""
    prefix, _, units = name.rpartition("_")
    return int(units) if prefix == "contig" and units.isdigit() else None


def secondary_mode_fraction(fit_metrics: dict[str, dict[str, Any]], min_gap: int) -> float | None:
    """Largest primary-record count >= ``min_gap`` units from the peak, relative to the peak."""
    counts: dict[int, int] = {}
    for contig, metrics in fit_metrics.items():
        units = _contig_units(contig)
        if units is not None:
            counts[units] = int(metrics.get("primary_alignment_records") or 0)
    if not counts or max(counts.values()) == 0:
        return None
    peak = max(counts, key=lambda contig: (counts[contig], -contig))
    distant = [n for contig, n in counts.items() if abs(contig - peak) >= min_gap]
    return max(distant, default=0) / counts[peak]


def assess_allele(info: dict[str, Any], settings: AlleleSelectionSettings) -> dict[str, Any]:
    """Return additive selection, depth and length-consistency fields for one allele."""
    fraction = secondary_mode_fraction(info.get("fit_metrics") or {}, settings.min_gap)
    unselected = (info.get("length_selection_evidence") or {}).get("unselected_passing_clusters")
    if fraction is not None and fraction >= settings.secondary_mode_min_fraction:
        selection = "unresolved_secondary_mode"
    elif unselected:
        selection = "unresolved_unselected_clusters"
    elif fraction is None:
        selection = "not_assessed"
    else:
        selection = "resolved"
    primary = info.get("primary_alignment_records")
    if primary is None:
        depth = "not_assessed"
    else:
        depth = "adequate" if primary >= settings.min_allele_primary_records else "low"
    length, reference_length = info.get("length"), info.get("reference_length")
    if length is None or reference_length is None:
        length_status = "not_assessed"
    elif length == reference_length:
        length_status = "consistent_with_consensus_contig"
    else:
        length_status = "cluster_center_differs_from_consensus_contig"
    return {
        "selection_status": selection,
        "secondary_mode_fraction": None if fraction is None else round(fraction, 4),
        "depth_status": depth,
        "depth_basis": "primary_alignment_records",
        "depth_threshold": settings.min_allele_primary_records,
        "length_status": length_status,
    }


def annotate_selection_qc(
    alleles: dict[str, Any], settings: AlleleSelectionSettings | None = None
) -> dict[str, Any]:
    """Add selection/depth/length gates to ``allele_1``/``allele_2`` in place."""
    settings = settings or DEFAULT_SETTINGS.allele_selection
    for key in ("allele_1", "allele_2"):
        info = alleles.get(key)
        if isinstance(info, dict):
            info.update(assess_allele(info, settings))
    return alleles
