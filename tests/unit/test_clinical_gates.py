"""Tests for the pure clinical evidence gates."""

from __future__ import annotations

from muc_one_span.clinical_gates import allele_gate_reasons, mutation_blockers

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
