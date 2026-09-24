"""Tests for the pure clinical evidence gates."""

from __future__ import annotations

from muc_one_span.clinical_gates import allele_gate_reasons, depth_gate_failure, mutation_blockers

_MISMATCH = "differs from the consensus contig length"


def test_length_gate_uses_recorded_length_status() -> None:
    """selection_qc's length_status is authoritative, not a recomputed comparison."""
    info = {"length": 44, "reference_length": 44}
    differs = {**info, "length_status": "cluster_center_differs_from_consensus_contig"}
    assert any(_MISMATCH in reason for reason in allele_gate_reasons(differs, "Allele 1"))
    consistent = {"length": 39, "reference_length": 44}
    consistent["length_status"] = "consistent_with_consensus_contig"
    assert not any(_MISMATCH in r for r in allele_gate_reasons(consistent, "Allele 1"))


def test_length_gate_falls_back_for_legacy_summaries() -> None:
    """Summaries without length_status keep the conservative length comparison."""
    reasons = allele_gate_reasons({"length": 39, "reference_length": 44}, "Allele 1")
    assert reasons == ["Allele 1: reported length 39 differs from the consensus contig length 44."]
    assert allele_gate_reasons({"length": 44, "reference_length": 44}, "Allele 1") == []


def test_unresolved_heterozygous_blocker_does_not_claim_an_allele_partition() -> None:
    """The status also arises on unphased same-length VCFs shared by both candidates."""
    mutation = {
        "frameshift": True,
        "template_match": True,
        "mutation_name": "dupC",
        "vcf_support": False,
        "vcf_support_status": "heterozygous_genotype_unresolved",
    }
    blockers = mutation_blockers(mutation)
    assert len(blockers) == 1
    assert "within its allele" not in blockers[0]
    assert blockers[0] == "heterozygous genotype not resolved to one allele"


_TEMPLATED = {
    "frameshift": True,
    "template_match": True,
    "mutation_name": "dupC",
    "localization_status": "exact",
    "vcf_support": False,
    "vcf_support_status": "not_applicable_read_consensus",
}


def test_read_support_statuses_name_the_blocker() -> None:
    for status in ("insufficient_depth", "discordant", "not_supported", "not_localized"):
        blockers = mutation_blockers({**_TEMPLATED, "read_support": {"status": status}})
        assert blockers == [f"read-level support {status}"]
    assert mutation_blockers({**_TEMPLATED, "read_support": {"status": "supported"}}) == []


def test_insufficient_depth_is_gated_like_low() -> None:
    info = {
        "depth_status": "insufficient",
        "depth_basis": "spanning_reads",
        "spanning_reads": 7,
        "depth_threshold": 30,
    }
    assert allele_gate_reasons(info, "Allele 1") == [
        "Allele 1: 7 spanning reads, below the per-allele depth gate (30)."
    ]


def test_depth_gate_fails_closed_on_unrecognised_status_with_a_basis() -> None:
    basis = {"depth_basis": "spanning_reads"}
    assert depth_gate_failure({**basis, "depth_status": "adequate"}) is None
    assert depth_gate_failure({**basis, "depth_status": "bogus"}, assessed=False) == "bogus"
    assert depth_gate_failure(basis) == "missing"
    assert depth_gate_failure({**basis, "depth_status": "not_assessed"}) == "not_assessed"
    assert depth_gate_failure({**basis, "depth_status": "not_assessed"}, assessed=False) is None
    assert depth_gate_failure({"depth_status": "bogus"}) is None  # legacy: no basis
    assert depth_gate_failure({"depth_status": "low"}, assessed=False) == "low"
    reasons = allele_gate_reasons({**basis, "depth_status": "bogus"}, "Allele 2")
    assert reasons == [
        "Allele 2: per-allele depth status 'bogus' (spanning reads) is not 'adequate'; "
        "the depth gate fails closed."
    ]
