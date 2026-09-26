"""Tests for clinical decision semantics under evidence and reconstruction limits."""

from __future__ import annotations

import pytest

from muc_one_span.report import compute_clinical_decision
from muc_one_span.settings import ClinicalDecisionSettings


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
            "template_match": True,
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


def test_heterozygous_length_partition_is_never_negative():
    """A length-partitioned allele with an unresolved heterozygous call is not NEGATIVE."""
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_2"].update(
        phase_status="single_heterozygous_unordered",
        allele_genotype_status="heterozygous_within_length_partition",
        consensus_haplotype="I",
        independent_haplotype_evidence=False,
    )
    assert compute_clinical_decision(summary)["state"] == "INCONCLUSIVE"


def _summary(mut: dict) -> dict:
    allele = {
        "length": 50,
        "canonical_repeats": 41,
        "reads": 200,
        "phase_status": "phased",
        "independent_haplotype_evidence": True,
    }
    return {
        "alleles": {
            "allele_1": dict(allele),
            "allele_2": dict(allele, length=60, canonical_repeats=51),
            "homozygous": False,
        },
        "classifications": {
            "allele_1": {
                "mutations": [mut],
                "ambiguous_bases": 0,
                "reconstruction_status": "complete_segmentation",
            },
            "allele_2": {
                "mutations": [],
                "ambiguous_bases": 0,
                "reconstruction_status": "complete_segmentation",
            },
        },
    }


BASE = {
    "repeat_index": 20,
    "mutation_name": "dupC",
    "frameshift": True,
    "localization_status": "exact",
    "template_match": True,
}


def test_missing_support_fields_are_not_pathogenic() -> None:
    decision = compute_clinical_decision(_summary(dict(BASE)))
    assert decision["state"] != "PATHOGENIC"


def test_read_support_is_accepted_as_explicit_support() -> None:
    mut = dict(
        BASE,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "supported", "alt": 180, "ref": 10, "other": 10},
    )
    assert compute_clinical_decision(_summary(mut))["state"] == "PATHOGENIC"


def test_vcf_exact_concordance_still_pathogenic() -> None:
    mut = dict(BASE, vcf_support=True, vcf_support_status="exact_sequence_concordance")
    assert compute_clinical_decision(_summary(mut))["state"] == "PATHOGENIC"


def test_unresolved_heterozygous_event_is_not_pathogenic() -> None:
    mut = dict(BASE, vcf_support=False, vcf_support_status="heterozygous_genotype_unresolved")
    assert compute_clinical_decision(_summary(mut))["state"] == "INCONCLUSIVE"


SUPPORTED = dict(BASE, vcf_support=True, vcf_support_status="exact_sequence_concordance")
RESOLVED_GATES = {
    "selection_status": "resolved",
    "secondary_mode_fraction": 0.01,
    "depth_status": "adequate",
    "depth_threshold": 30,
    "primary_alignment_records": 400,
    "length_status": "consistent_with_consensus_contig",
}


def _gated_summary(mutations: list[dict] | None = None) -> dict:
    summary = _summary(dict(BASE))
    summary["classifications"]["allele_1"]["mutations"] = mutations or []
    summary["run_status"] = {"status": "completed"}
    for key in ("allele_1", "allele_2"):
        summary["alleles"][key].update(RESOLVED_GATES)
    return summary


def test_resolved_gates_keep_negative() -> None:
    decision = compute_clinical_decision(_gated_summary())
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"


def test_secondary_mode_blocks_negative() -> None:
    """MP3-shaped: a second primary-record mode at 0.62 of the peak inside one cluster."""
    summary = _gated_summary()
    summary["alleles"]["allele_1"].update(
        selection_status="unresolved_secondary_mode", secondary_mode_fraction=0.619
    )
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any("allele selection unresolved" in detail for detail in decision["details"])


def test_low_allele_depth_blocks_negative_and_pathogenic() -> None:
    negative = _gated_summary()
    negative["alleles"]["allele_2"].update(depth_status="low", primary_alignment_records=12)
    assert compute_clinical_decision(negative)["state"] == "INCONCLUSIVE"
    positive = _gated_summary([dict(SUPPORTED)])
    positive["alleles"]["allele_1"].update(depth_status="low", primary_alignment_records=12)
    decision = compute_clinical_decision(positive)
    assert decision["state"] == "INCONCLUSIVE"
    assert any("per-allele depth gate" in detail for detail in decision["details"])


def test_legacy_low_total_reads_blocks_pathogenic() -> None:
    summary = _summary(dict(SUPPORTED))
    summary["alleles"]["allele_1"]["reads"] = 10
    summary["alleles"]["allele_2"]["reads"] = 12
    assert compute_clinical_decision(summary)["state"] == "INCONCLUSIVE"


def test_untemplated_frameshift_is_not_pathogenic() -> None:
    """MP5-shaped: a supported novel frameshift without a dictionary template is inconclusive."""
    novel = {
        "repeat_index": 35,
        "closest_type": "X",
        "differences": [{"type": "insertion", "position": 23}],
        "frameshift": True,
        "vcf_support": True,
        "vcf_support_status": "exact_sequence_concordance",
        "localization_status": "resolved",
    }
    decision = compute_clinical_decision(_gated_summary([novel]))
    assert decision["state"] == "INCONCLUSIVE"
    assert any("event identity not established" in detail for detail in decision["details"])


def test_length_contig_mismatch_blocks_negative() -> None:
    """MP1-shaped allele 1: reported 39 units, consensus from the 44-unit contig."""
    summary = _gated_summary()
    summary["alleles"]["allele_1"].update(
        length=39,
        reference_length=44,
        length_status="cluster_center_differs_from_consensus_contig",
    )
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any(
        "reported length 39 differs from the consensus contig length 44" in detail
        for detail in decision["details"]
    )


def test_pathogenic_lists_quality_caveats() -> None:
    summary = _gated_summary([dict(SUPPORTED)])
    summary["alleles"]["allele_2"].update(
        selection_status="unresolved_secondary_mode", secondary_mode_fraction=0.28
    )
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "PATHOGENIC"
    assert decision["details"][0].startswith("Allele 1: dupC at repeat unit 20")
    assert any(detail.startswith("Quality caveat: Allele 2") for detail in decision["details"])


def test_one_low_allele_without_depth_on_other_blocks_negative() -> None:
    summary = _gated_summary()
    for key in ("depth_status", "depth_threshold", "primary_alignment_records"):
        summary["alleles"]["allele_1"].pop(key)
    summary["alleles"]["allele_2"].update(depth_status="low", primary_alignment_records=12)
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any("Allele 2: 12 primary alignments" in detail for detail in decision["details"])


def test_not_assessed_depth_falls_back_to_legacy_total_reads() -> None:
    negative = _gated_summary()
    positive = _gated_summary([dict(SUPPORTED)])
    for summary in (negative, positive):
        for key in ("allele_1", "allele_2"):
            summary["alleles"][key].update(depth_status="not_assessed", reads=10)
    decision = compute_clinical_decision(negative)
    assert decision["state"] == "INCONCLUSIVE"
    assert any("Total read depth (20 reads)" in detail for detail in decision["details"])
    assert compute_clinical_decision(positive)["state"] == "INCONCLUSIVE"


def test_pathogenic_allele_1_with_low_allele_2_keeps_caveat() -> None:
    summary = _gated_summary([dict(SUPPORTED)])
    summary["alleles"]["allele_2"].update(depth_status="low", primary_alignment_records=12)
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "PATHOGENIC"
    assert any(
        detail.startswith("Quality caveat: Allele 2: 12 primary alignments")
        for detail in decision["details"]
    )


HYBRID_GATES = dict(
    RESOLVED_GATES,
    depth_basis="spanning_reads",
    spanning_reads=120,
    secondary_mode_fraction=None,
    allele_genotype_status="not_applicable_read_consensus",
    engine="hybrid",
)


def _hybrid_summary(**allele2: object) -> dict:
    summary = _gated_summary()
    for key in ("allele_1", "allele_2"):
        summary["alleles"][key].update(HYBRID_GATES)
    summary["alleles"]["allele_2"].update(allele2)
    return summary


def test_hybrid_resolved_is_negative() -> None:
    assert compute_clinical_decision(_hybrid_summary())["state"] == (
        "NO_PATHOGENIC_VARIANT_DETECTED"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("depth_status", "insufficient"),
        ("selection_status", "unresolved_rejected_peak"),
        ("selection_status", "unresolved_unassigned_spanning"),
        ("selection_status", "unresolved_single_site"),
        ("selection_status", "unresolved_max_alleles"),
        ("allele_genotype_status", "residual_heterogeneity"),
    ],
)
def test_hybrid_gate_blocks_negative(field: str, value: str) -> None:
    decision = compute_clinical_decision(_hybrid_summary(**{field: value}))
    assert decision["state"] == "INCONCLUSIVE"


def test_insufficient_depth_carrier_blocks_pathogenic() -> None:
    mutation = dict(
        BASE,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "supported"},
    )
    summary = _hybrid_summary()
    summary["classifications"]["allele_1"]["mutations"] = [mutation]
    assert compute_clinical_decision(summary)["state"] == "PATHOGENIC"
    summary["alleles"]["allele_1"]["depth_status"] = "insufficient"
    assert compute_clinical_decision(summary)["state"] == "INCONCLUSIVE"


def test_hybrid_depth_message_names_spanning_reads() -> None:
    decision = compute_clinical_decision(_hybrid_summary(depth_status="low", spanning_reads=12))
    assert any("Allele 2: 12 spanning reads" in detail for detail in decision["details"])


def test_selection_message_omits_missing_secondary_fraction() -> None:
    decision = compute_clinical_decision(
        _hybrid_summary(selection_status="unresolved_rejected_peak", selection_detail="peak 80u")
    )
    reason = next(d for d in decision["details"] if "allele selection unresolved" in d)
    assert "secondary mode fraction" not in reason and "peak 80u" in reason


def test_unknown_depth_status_fails_closed_when_other_allele_is_adequate() -> None:
    decision = compute_clinical_decision(_hybrid_summary(depth_status="bogus"))
    assert decision["state"] == "INCONCLUSIVE"
    assert any("Allele 2" in d and "'bogus'" in d for d in decision["details"])


def test_not_assessed_depth_beside_an_assessed_allele_fails_closed() -> None:
    negative = _hybrid_summary(depth_status="not_assessed")
    assert compute_clinical_decision(negative)["state"] == "INCONCLUSIVE"
    for status in ("not_assessed", "bogus"):
        positive = _hybrid_summary()
        positive["classifications"]["allele_1"]["mutations"] = [dict(SUPPORTED)]
        positive["alleles"]["allele_1"]["depth_status"] = status
        decision = compute_clinical_decision(positive)
        assert decision["state"] == "INCONCLUSIVE"
        expected = f"carrying allele depth status {status!r} is not adequate"
        assert any(expected in d for d in decision["details"])


def test_ladder_not_assessed_on_both_alleles_keeps_legacy_fallback() -> None:
    """Ladder without a BAM marks both alleles not_assessed; total reads decide."""
    summary = _gated_summary()
    for key in ("allele_1", "allele_2"):
        summary["alleles"][key].update(
            depth_status="not_assessed", depth_basis="primary_alignment_records", reads=40
        )
    assert compute_clinical_decision(summary)["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"


def test_non_string_depth_basis_fails_closed() -> None:
    summary = _hybrid_summary(depth_basis=["spanning_reads"])
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any("Allele 2" in d and "invalid_depth_basis" in d for d in decision["details"])


_R9_DISCORDANCE = {
    "status": "discordant_frameshift",
    "min_af": 0.5,
    "min_depth": 10,
    "records": [
        {"chrom": "contig_82", "pos": 2833, "ref": "C", "alt": "CG", "af": 0.719, "dp": 121}
    ],
}


def test_stage_discordance_blocks_negative() -> None:
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_2"]["stage_concordance"] = dict(_R9_DISCORDANCE)
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any(d.startswith("Allele 2: caller-stage discordance:") for d in decision["details"])


def test_stage_discordance_is_caveat_on_pathogenic() -> None:
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_2"]["stage_concordance"] = dict(_R9_DISCORDANCE)
    summary["classifications"]["allele_1"]["mutations"] = [
        {
            "mutation_name": "dupC",
            "repeat_index": 17,
            "frameshift": True,
            "template_match": True,
            "vcf_support": True,
            "vcf_support_status": "exact_sequence_concordance",
            "localization_status": "resolved",
        }
    ]
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "PATHOGENIC"
    assert any(
        d.startswith("Quality caveat: Allele 2: caller-stage discordance:")
        for d in decision["details"]
    )


def test_concordant_or_missing_stage_record_keeps_v0160_negative() -> None:
    baseline = compute_clinical_decision(_resolved_diploid_summary())
    assert baseline["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"
    summary = _resolved_diploid_summary()
    for key in ("allele_1", "allele_2"):
        summary["alleles"][key]["stage_concordance"] = {"status": "concordant", "records": []}
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"
    assert decision["details"] == baseline["details"]


def test_not_assessed_stage_record_keeps_negative_with_caveat() -> None:
    """A missing pileup VCF is not silent (R16) but does not block NEGATIVE."""
    baseline = compute_clinical_decision(_resolved_diploid_summary())
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_1"]["stage_concordance"] = {
        "status": "not_assessed",
        "reason": "pileup_vcf_unavailable",
        "records": [],
    }
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"
    assert decision["details"][: len(baseline["details"])] == baseline["details"]
    assert decision["details"][len(baseline["details"]) :] == [
        "Quality caveat: Allele 1: caller-stage concordance not assessed "
        "(pileup_vcf_unavailable); Clair3 pileup-stage frameshift calls were not compared "
        "with the applied calls."
    ]


def test_not_assessed_stage_record_is_caveat_on_inconclusive() -> None:
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_1"]["stage_concordance"] = {
        "status": "not_assessed",
        "reason": "pileup_vcf_unavailable",
    }
    summary["alleles"]["allele_2"]["stage_concordance"] = dict(_R9_DISCORDANCE)
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert decision["details"][-1].startswith(
        "Quality caveat: Allele 1: caller-stage concordance not assessed"
    )


def test_ambiguous_bases_at_default_threshold_stays_negative() -> None:
    """Summed ambiguous bases at the v0.16.0 default (10) must not block NEGATIVE."""
    summary = _gated_summary()
    summary["classifications"]["allele_1"]["ambiguous_bases"] = 6
    summary["classifications"]["allele_2"]["ambiguous_bases"] = 4
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"
    assert decision["thresholds"] == {
        "max_ambiguous_bases": 10,
        "legacy_min_total_reads": 30,
        "source": "default",
    }


def test_ambiguous_bases_above_default_threshold_is_inconclusive() -> None:
    """Byte-identical to the v0.16.0 wording at the default threshold."""
    summary = _gated_summary()
    summary["classifications"]["allele_1"]["ambiguous_bases"] = 6
    summary["classifications"]["allele_2"]["ambiguous_bases"] = 5
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert "High number of ambiguous consensus bases (11) detected." in decision["details"]


def test_recorded_max_ambiguous_bases_raises_the_threshold() -> None:
    summary = _gated_summary()
    summary["classifications"]["allele_1"]["ambiguous_bases"] = 15
    summary["configuration"] = {
        "settings": {"clinical_decision": {"max_ambiguous_bases": 20, "legacy_min_total_reads": 30}}
    }
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"
    assert decision["thresholds"]["source"] == "recorded_configuration"


def test_explicit_settings_override_recorded_configuration() -> None:
    summary = _gated_summary()
    summary["classifications"]["allele_1"]["ambiguous_bases"] = 15
    summary["configuration"] = {
        "settings": {"clinical_decision": {"max_ambiguous_bases": 20, "legacy_min_total_reads": 30}}
    }
    decision = compute_clinical_decision(
        summary, settings=ClinicalDecisionSettings(max_ambiguous_bases=5)
    )
    assert decision["state"] == "INCONCLUSIVE"
    assert decision["thresholds"]["source"] == "explicit"


def test_legacy_low_total_reads_detail_text_is_byte_identical() -> None:
    summary = _summary(dict(SUPPORTED))
    summary["alleles"]["allele_1"]["reads"] = 20
    summary["alleles"]["allele_2"]["reads"] = 5
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert (
        "Total read depth (25 reads) is below diagnostic threshold (30 reads)."
        in decision["details"]
    )


def test_recorded_legacy_min_total_reads_avoids_low_coverage() -> None:
    summary = _summary(dict(SUPPORTED))
    summary["alleles"]["allele_1"]["reads"] = 20
    summary["alleles"]["allele_2"]["reads"] = 5
    summary["configuration"] = {"settings": {"clinical_decision": {"legacy_min_total_reads": 20}}}
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "PATHOGENIC"
    assert not any("below diagnostic threshold" in detail for detail in decision["details"])


def test_stage_veto_and_hybrid_gates_both_apply() -> None:
    """Merged decision: the stage veto and the hybrid gates each block NEGATIVE only."""
    vetoed = _hybrid_summary()
    vetoed["alleles"]["allele_2"]["stage_concordance"] = dict(_R9_DISCORDANCE)
    decision = compute_clinical_decision(vetoed)
    assert decision["state"] == "INCONCLUSIVE"
    assert any(d.startswith("Allele 2: caller-stage discordance:") for d in decision["details"])
    both = _hybrid_summary(selection_status="unresolved_single_site")
    both["alleles"]["allele_2"]["stage_concordance"] = dict(_R9_DISCORDANCE)
    assert compute_clinical_decision(both)["state"] == "INCONCLUSIVE"
    carrier = _hybrid_summary()
    carrier["alleles"]["allele_2"]["stage_concordance"] = dict(_R9_DISCORDANCE)
    carrier["classifications"]["allele_1"]["mutations"] = [dict(SUPPORTED)]
    assert compute_clinical_decision(carrier)["state"] == "PATHOGENIC"
