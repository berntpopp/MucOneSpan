"""Audit: every hybrid status that leaves an allele unaccounted blocks NEGATIVE.

Each summary here is built with the engine's own record builders (``allele_info``,
``single_group_fields``), so the audit follows the producer, and every decision goes
through ``compute_clinical_decision``. The benchmark's ``missing_allele``,
``unresolved_allele_alias`` and ``ambiguous_reconstruction`` are evaluator labels
(truth comparison, and the evaluator's reading of every one-group allele_2 alias); the
caller never emits them, so they cannot be gates. What the caller emits for the same
situations is audited below.
"""

from __future__ import annotations

from typing import Any

import pytest

from muc_one_span.hybrid.allele_fields import (
    PHASE_STATUS,
    RESOLVED,
    UNRESOLVED_SELECTION_STATUSES,
    allele_info,
    single_group_fields,
)
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.report import compute_clinical_decision
from muc_one_span.settings import DEFAULT_SETTINGS

H = DEFAULT_SETTINGS.hybrid
UNIT = 60
FIXED = DEFAULT_SETTINGS.reference_layout.fixed_repeat_count
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
RESIDUAL = [{"pos": 100, "af": 0.3, "kind": "col"}]


def _members(n: int) -> list[SpanRead]:
    return [SpanRead(f"r{i}", "A", 20.0, "+", 0, "motif") for i in range(n)]


def _allele(
    name: str,
    units: int,
    *,
    basis: str = "length",
    selection: str = RESOLVED,
    reads: int = H.depth_adequate_spanning,
    residual: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    return allele_info(
        name,
        "A" * units * UNIT,
        _members(reads),
        0,
        basis,
        residual or [],
        (selection, "detail"),
        UNIT,
        FIXED,
        H,
        concordance=1.0,
    )


def _summary(alleles: dict[str, Any], mutations: list[dict[str, Any]] | None = None) -> dict:
    classes = {
        k: {"mutations": [], "ambiguous_bases": 0, "reconstruction_status": "complete_segmentation"}
        for k in alleles
        if isinstance(alleles[k], dict) and "candidate_duplicate_of" not in alleles[k]
    }
    if mutations:
        classes["allele_1"]["mutations"] = mutations
    return {"run_status": {"status": "completed"}, "alleles": alleles, "classifications": classes}


def _two(**kw: Any) -> dict[str, Any]:
    alleles = {"allele_1": _allele("allele_1", 40, **kw), "allele_2": _allele("allele_2", 60, **kw)}
    alleles.update(homozygous=False, same_length=False, allele_multiplicity_status="resolved")
    return alleles


def _one(**kw: Any) -> dict[str, Any]:
    """One allele group, as the engine reports it (``single_group_fields``)."""
    kw.setdefault("basis", "none")
    first = _allele("allele_1", 60, **kw)
    residual_any = bool(first["residual_sites"])
    return single_group_fields(first, first["selection_status"], residual_any)


def test_controls_are_negative() -> None:
    """Resolved two-length and resolved homozygous summaries stay NEGATIVE."""
    assert compute_clinical_decision(_summary(_two()))["state"] == NEGATIVE
    homozygous = _one()
    assert homozygous["homozygous"] is True
    assert compute_clinical_decision(_summary(homozygous))["state"] == NEGATIVE


@pytest.mark.parametrize("status", sorted(UNRESOLVED_SELECTION_STATUSES))
@pytest.mark.parametrize("groups", ["one", "two"])
def test_every_unresolved_selection_status_blocks_negative(status: str, groups: str) -> None:
    alleles = _one(selection=status) if groups == "one" else _two(selection=status)
    decision = compute_clinical_decision(_summary(alleles))
    assert decision["state"] == "INCONCLUSIVE", status
    assert any(status in d for d in decision["details"]), decision["details"]


def test_selection_statuses_cover_every_unconfirmed_phase_status() -> None:
    unresolved = {v for v in PHASE_STATUS.values() if v.startswith("unresolved")}
    assert unresolved <= UNRESOLVED_SELECTION_STATUSES


@pytest.mark.parametrize("basis", sorted(b for b, s in PHASE_STATUS.items() if "unresolved" in s))
def test_unconfirmed_one_group_is_not_homozygous(basis: str) -> None:
    """An unconfirmed phase basis never yields a homozygous (resolved) one-group call."""
    status = PHASE_STATUS[basis]
    alleles = _one(basis=basis, selection=status)
    assert alleles["homozygous"] is False
    assert alleles["allele_multiplicity_status"] == "unresolved"
    assert alleles["sequence_identity_status"] == "unresolved"
    assert alleles["allele_2"]["independent_haplotype_evidence"] is False
    assert compute_clinical_decision(_summary(alleles))["state"] == "INCONCLUSIVE"


def test_one_group_with_residual_heterogeneity_is_inconclusive() -> None:
    alleles = _one(residual=RESIDUAL)
    assert alleles["homozygous"] is False
    decision = compute_clinical_decision(_summary(alleles))
    assert decision["state"] == "INCONCLUSIVE"
    assert any("residual read heterogeneity" in d for d in decision["details"])


def test_two_groups_with_residual_heterogeneity_is_inconclusive() -> None:
    assert compute_clinical_decision(_summary(_two(residual=RESIDUAL)))["state"] == "INCONCLUSIVE"


@pytest.mark.parametrize("reads", [H.depth_low_spanning, H.depth_low_spanning - 1])
def test_low_or_insufficient_allele_depth_is_inconclusive(reads: int) -> None:
    alleles = _two()
    alleles["allele_2"] = _allele("allele_2", 60, reads=reads)
    assert alleles["allele_2"]["depth_status"] != "adequate"
    assert compute_clinical_decision(_summary(alleles))["state"] == "INCONCLUSIVE"


def test_unknown_depth_status_fails_closed() -> None:
    alleles = _two()
    alleles["allele_2"]["depth_status"] = "unknown"
    assert compute_clinical_decision(_summary(alleles))["state"] == "INCONCLUSIVE"


@pytest.mark.parametrize(
    "status", ["insufficient_depth", "discordant", "not_supported", "not_localized"]
)
def test_unsupported_event_is_never_negative_or_pathogenic(status: str) -> None:
    mutation = {
        "repeat_index": 20,
        "mutation_name": "dupC",
        "frameshift": True,
        "localization_status": "exact",
        "template_match": True,
        "vcf_support": False,
        "vcf_support_status": "not_applicable_read_consensus",
        "read_support": {"status": status},
    }
    decision = compute_clinical_decision(_summary(_two(), [mutation]))
    assert decision["state"] == "INCONCLUSIVE", status
