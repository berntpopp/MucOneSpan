"""Scientific regression contracts for classification completeness and frame."""

from muc_one_span.classify import classify_repeat, classify_sequence, validate_mutations_against_vcf
from muc_one_span.config import load_repeat_dictionary


def test_signed_net_frame_for_mixed_indels() -> None:
    rd = load_repeat_dictionary()
    x = rd.repeats["X"]
    result = classify_repeat(x[:5] + x[7:30] + "A" + x[30:], rd)
    assert result["frameshift"] is True
    assert result["net_indel_bases"] == -1


def test_balanced_indels_restore_downstream_frame() -> None:
    rd = load_repeat_dictionary()
    x = rd.repeats["X"]
    result = classify_repeat(x[:5] + x[6:30] + "A" + x[30:], rd)
    assert result["frameshift"] is False
    assert result["net_indel_bases"] == 0


def test_unclassified_tail_is_in_denominator() -> None:
    rd = load_repeat_dictionary()
    result = classify_sequence(rd.repeats["X"] * 5 + "N" * 29, rd)
    assert result["unclassified_regions"] == [
        {"start": 300, "end": 329, "reason": "unresolved_sequence"}
    ]
    assert result["classified_bases"] == 300
    assert result["sequence_length"] == 329
    assert result["ambiguous_bases"] == 29
    assert result["reconstruction_status"] == "ambiguous_reconstruction"
    assert result["classification_coverage"] == 300 / 329


def test_large_unresolved_insertion_does_not_invent_mutation_or_repeat_index() -> None:
    rd = load_repeat_dictionary()
    x = rd.repeats["X"]
    seq = x * 2 + "A" * 40 + x * 3 + rd.repeats["9"]
    result = classify_sequence(seq, rd, strict_segmentation=True)
    assert result["reconstruction_status"] == "ambiguous_reconstruction"
    assert result["mutations_detected"] == []
    spans = [(r["start"], r["end"]) for r in result["repeats"]]
    spans += [(r["start"], r["end"]) for r in result["unclassified_regions"]]
    spans += [(r["start"], r["end"]) for r in result["recovered_suffix"]]
    assert sorted(spans)[0][0] == 0
    assert sorted(spans)[-1][1] == len(seq)
    assert all(a[1] == b[0] for a, b in zip(sorted(spans), sorted(spans)[1:], strict=False))
    assert all(r["index"] is None for r in result["recovered_suffix"])


def test_position_only_record_is_not_variant_identity_support() -> None:
    rd = load_repeat_dictionary()
    mutant = next(s for s, label in rd.mutated_sequences.items() if label == ("X", "dupC"))
    result = classify_sequence(rd.repeats["X"] + mutant, rd)
    validated = validate_mutations_against_vcf(result, [{"pos": 585, "qual": 40}])
    assert validated["mutations_detected"][0]["vcf_support"] is False
    assert validated["mutations_detected"][0]["vcf_support_status"] == "projection_unavailable"


def test_candidate_default_retains_downstream_templates_after_noisy_window() -> None:
    rd = load_repeat_dictionary()
    x = rd.repeats["X"]
    mutant = next(s for s, label in rd.mutated_sequences.items() if label == ("X", "dupC"))
    seq = x + "N" * 60 + x + mutant + x
    result = classify_sequence(seq, rd)
    assert any(m.get("mutation_name") == "dupC" for m in result["mutations_detected"])
    assert result["reconstruction_status"] == "ambiguous_reconstruction"
    assert all(m["localization_status"] == "ambiguous" for m in result["mutations_detected"])


def test_projection_unavailable_preserves_fit_confidence_while_absent_penalizes(
    monkeypatch,
) -> None:
    classification = {
        "mutations_detected": [{"repeat_index": 1}],
        "repeats": [{"confidence": 0.8}],
    }

    monkeypatch.setattr(
        "muc_one_span.variant_support.mutation_concordance",
        lambda *_: {1: ("projection_unavailable", [])},
    )
    unavailable = validate_mutations_against_vcf(classification, [])
    mutation = unavailable["mutations_detected"][0]
    assert mutation["vcf_support"] is False
    assert mutation["vcf_support_status"] == "projection_unavailable"
    assert mutation["vcf_qual"] == 0.0
    assert unavailable["repeats"][0]["confidence"] == 0.8
    assert unavailable["allele_confidence"] == 0.8

    monkeypatch.setattr(
        "muc_one_span.variant_support.mutation_concordance",
        lambda *_: {1: ("localization_ambiguous", [])},
    )
    ambiguous = validate_mutations_against_vcf(classification, [])
    assert ambiguous["mutations_detected"][0]["vcf_support_status"] == "localization_ambiguous"
    assert ambiguous["repeats"][0]["confidence"] == 0.8
    assert ambiguous["allele_confidence"] == 0.8

    monkeypatch.setattr(
        "muc_one_span.variant_support.mutation_concordance",
        lambda *_: {1: ("absent", [])},
    )
    absent = validate_mutations_against_vcf(classification, [])
    assert absent["mutations_detected"][0]["vcf_support_status"] == "absent"
    assert absent["repeats"][0]["confidence"] == 0.24
    assert absent["allele_confidence"] == 0.24
