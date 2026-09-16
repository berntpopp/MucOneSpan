"""Tests for clinical decision semantics under evidence and reconstruction limits."""

from __future__ import annotations

from muc_one_span.report import compute_clinical_decision


def _resolved_diploid_summary():
    return {
        "alleles": {
            "allele_1": {
                "length": 50,
                "reads": 500,
                "canonical_repeats": 41,
                "contig_name": "contig_41",
                "independent_haplotype_evidence": True,
                "phase_status": "phased",
            },
            "allele_2": {
                "length": 65,
                "reads": 400,
                "canonical_repeats": 56,
                "contig_name": "contig_56",
                "independent_haplotype_evidence": True,
                "phase_status": "phased",
            },
            "homozygous": False,
        },
        "classifications": {
            "allele_1": {"mutations": [], "ambiguous_bases": 0},
            "allele_2": {"mutations": [], "ambiguous_bases": 0},
        },
        "run_status": {"status": "completed"},
    }


def test_supported_frameshift_yields_pathogenic():
    summary = _resolved_diploid_summary()
    summary["classifications"]["allele_1"]["mutations"] = [
        {
            "mutation_name": "dupC",
            "repeat_index": 17,
            "frameshift": True,
            "vcf_support": True,
            "vcf_support_status": "exact_sequence_concordance",
            "localization_status": "resolved",
        }
    ]
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "PATHOGENIC"
    assert "Allele 1: dupC at repeat unit 17" in decision["details"][0]


def test_in_frame_ambiguous_expansion_yields_inconclusive():
    """HG002 WGS-style in-frame 18bp expansion must not produce PATHOGENIC."""
    summary = _resolved_diploid_summary()
    summary["classifications"]["allele_2"]["mutations"] = [
        {
            "mutation_name": "ins18",
            "repeat_index": 70,
            "frameshift": False,
            "vcf_support": False,
            "vcf_support_status": "absent",
            "localization_status": "ambiguous",
        }
    ]
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any(
        "inconclusive" in d.lower() or "unresolved" in d.lower() or "ins18" in d
        for d in decision["details"]
    )


def test_unsupported_frameshift_mutation_yields_inconclusive():
    """Frameshift with absent VCF support must not produce definitive PATHOGENIC."""
    summary = _resolved_diploid_summary()
    summary["classifications"]["allele_1"]["mutations"] = [
        {
            "mutation_name": "dupC",
            "repeat_index": 80,
            "frameshift": True,
            "vcf_support": False,
            "vcf_support_status": "absent",
            "localization_status": "ambiguous",
        }
    ]
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"


def test_unresolved_candidate_reconstruction_yields_inconclusive():
    """HG002 PCR-style unresolved second candidate must not produce a reassuring NEGATIVE."""
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_2"] = {
        "length": 15,
        "reads": 100,
        "canonical_repeats": 6,
        "contig_name": "contig_6",
        "candidate_duplicate_of": "allele_1",
        "independent_haplotype_evidence": False,
        "reconstruction_status": "candidate_reference_confidence_unverified",
        "phase_status": "unphased",
    }
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any(
        "reconstruction" in d.lower() or "unresolved" in d.lower() or "phase" in d.lower()
        for d in decision["details"]
    )


def test_unphased_candidate_yields_inconclusive():
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_1"]["independent_haplotype_evidence"] = False
    summary["alleles"]["allele_1"]["phase_status"] = "unphased"
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"


def test_resolved_negative_yields_no_pathogenic_variant_detected():
    summary = _resolved_diploid_summary()
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"
