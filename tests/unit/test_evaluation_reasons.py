"""Clinical reason lists and evaluator reconstruction flags kept per evaluated sample."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from muc_one_span.evaluation.clinical_confusion import predicted_clinical, predicted_decision
from muc_one_span.evaluation.models import PredictedAllele, RunObservation
from muc_one_span.evaluation.reasons import reconstruction_flags
from muc_one_span.report import compute_clinical_decision
from muc_one_span.settings import DEFAULT_SETTINGS

GATE = DEFAULT_SETTINGS.allele_selection.min_allele_primary_records


def _allele(length: int, records: int) -> dict[str, Any]:
    status = "adequate" if records >= GATE else "low"
    return {
        "length": length,
        "reads": records,
        "primary_alignment_records": records,
        "depth_status": status,
        "depth_threshold": GATE,
    }


def _summary(records: int) -> dict[str, Any]:
    return {
        "alleles": {"allele_1": _allele(40, records), "allele_2": _allele(41, GATE)},
        "classifications": {},
    }


def test_predicted_clinical_keeps_the_inconclusive_reason_list(tmp_path: Path) -> None:
    summary = _summary(GATE - 1)
    (tmp_path / "summary.json").write_text(json.dumps(summary))
    result = predicted_clinical(tmp_path)
    assert result["decision"] == "INCONCLUSIVE" == predicted_decision(tmp_path)
    assert result["reasons"] == compute_clinical_decision(summary)["details"]
    assert any("per-allele depth gate" in r for r in result["reasons"])


def test_predicted_clinical_negative_and_unreadable_have_no_reasons(tmp_path: Path) -> None:
    assert predicted_clinical(tmp_path) == {"decision": "NO_CALL", "reasons": []}
    (tmp_path / "summary.json").write_text(json.dumps(_summary(GATE)))
    assert predicted_clinical(tmp_path) == {
        "decision": "NO_PATHOGENIC_VARIANT_DETECTED",
        "reasons": [],
    }


def _pred(name: str, sequence: str) -> PredictedAllele:
    return PredictedAllele(name, sequence, ("X",), 1, 1)


def test_reconstruction_flags_name_each_evaluator_cause() -> None:
    iupac = RunObservation("ambiguous_reconstruction", (_pred("a", "ACGR"),))
    row = {"missing_alleles": 1, "unproven_duplicate_alleles": 0, "extra_alleles": 0}
    assert reconstruction_flags(iupac, row) == [
        "ambiguous_reconstruction",
        "iupac_bases",
        "missing_allele",
    ]
    alias = RunObservation(
        "ambiguous_reconstruction",
        (_pred("a", "ACGT"),),
        warnings=("unresolved_allele_alias:allele_2:allele_1",),
    )
    assert reconstruction_flags(alias, {"unproven_duplicate_alleles": 1}) == [
        "ambiguous_reconstruction",
        "unresolved_allele_alias",
        "unproven_duplicate_allele",
    ]
    status = RunObservation("ambiguous_reconstruction", (_pred("a", "ACGT"),))
    assert reconstruction_flags(status, {"extra_alleles": 1}) == [
        "ambiguous_reconstruction",
        "producer_status_unresolved",
        "extra_allele",
    ]
    assert reconstruction_flags(RunObservation("insufficient_evidence"), {}) == [
        "insufficient_evidence"
    ]
    assert reconstruction_flags(RunObservation("completed", (_pred("a", "ACGT"),)), {}) == []
