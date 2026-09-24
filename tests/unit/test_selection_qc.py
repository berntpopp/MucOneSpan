"""Allele-selection, depth and length gates from ladder primary-alignment evidence."""

from __future__ import annotations

import pytest

from muc_one_span.selection_qc import annotate_selection_qc, assess_allele, secondary_mode_fraction
from muc_one_span.settings import AlleleSelectionSettings


def _metrics(counts: dict[int, int]) -> dict[str, dict[str, int]]:
    return {f"contig_{unit}": {"primary_alignment_records": n} for unit, n in counts.items()}


# Shapes of PRJEB92208 cohort-v3 clusters (public libraries): a control peak with a
# 1% distant tail, and MP3 allele 1 merging a second allele mode at 0.62 of the peak.
CONTROL = _metrics({6: 36, 27: 284, 28: 2855, 29: 3510, 30: 392})
MERGED = _metrics({32: 2754, 33: 2872, 34: 400, 38: 1778, 39: 1679})


def test_secondary_mode_fraction_separates_control_and_merged_cluster() -> None:
    assert secondary_mode_fraction(CONTROL, 5) == pytest.approx(36 / 3510)
    assert secondary_mode_fraction(MERGED, 5) == pytest.approx(1778 / 2872)
    assert secondary_mode_fraction({}, 5) is None
    assert secondary_mode_fraction(_metrics({10: 0}), 5) is None


def test_assess_allele_reports_all_gates() -> None:
    settings = AlleleSelectionSettings()
    base = {"length": 45, "reference_length": 45, "primary_alignment_records": 632}
    assert assess_allele({**base, "fit_metrics": CONTROL}, settings) == {
        "selection_status": "resolved",
        "secondary_mode_fraction": round(36 / 3510, 4),
        "depth_status": "adequate",
        "depth_basis": "primary_alignment_records",
        "depth_threshold": 30,
        "length_status": "consistent_with_consensus_contig",
    }
    merged = assess_allele({**base, "fit_metrics": MERGED}, settings)
    assert merged["selection_status"] == "unresolved_secondary_mode"
    low = assess_allele({**base, "fit_metrics": CONTROL, "primary_alignment_records": 12}, settings)
    assert low["depth_status"] == "low"
    shifted = assess_allele(
        {**base, "length": 39, "reference_length": 44, "fit_metrics": CONTROL}, settings
    )
    assert shifted["length_status"] == "cluster_center_differs_from_consensus_contig"
    unassessed = assess_allele({"length": 45}, settings)
    assert (
        unassessed["selection_status"],
        unassessed["depth_status"],
        unassessed["length_status"],
    ) == ("not_assessed", "not_assessed", "not_assessed")


def test_unselected_passing_cluster_marks_selection_unresolved() -> None:
    info = {
        "fit_metrics": CONTROL,
        "primary_alignment_records": 632,
        "length_selection_evidence": {"unselected_passing_clusters": 1},
    }
    status = assess_allele(info, AlleleSelectionSettings())["selection_status"]
    assert status == "unresolved_unselected_clusters"


def test_gate_thresholds_come_from_settings() -> None:
    settings = AlleleSelectionSettings(
        secondary_mode_min_fraction=0.7, min_allele_primary_records=700
    )
    result = assess_allele({"fit_metrics": MERGED, "primary_alignment_records": 632}, settings)
    assert (result["selection_status"], result["depth_status"]) == ("resolved", "low")
    assert result["depth_threshold"] == 700


def test_annotate_updates_both_alleles_in_place() -> None:
    alleles = {
        "allele_1": {"fit_metrics": CONTROL, "primary_alignment_records": 632},
        "allele_2": {"fit_metrics": MERGED, "primary_alignment_records": 80},
        "homozygous": False,
    }
    assert annotate_selection_qc(alleles) is alleles
    assert alleles["allele_1"]["selection_status"] == "resolved"
    assert alleles["allele_2"]["selection_status"] == "unresolved_secondary_mode"
    assert alleles["homozygous"] is False
