"""Adversarial contracts for strict offline scoring."""

from muc_one_span.evaluation.matching import match_alleles
from muc_one_span.evaluation.models import (
    Event,
    PredictedAllele,
    RunObservation,
    TruthHaplotype,
    TruthSample,
)
from muc_one_span.evaluation.scoring import aggregate, evaluate_sample, score_events


def hap(name="h1", seq="A", events=()):
    return TruthHaplotype(name, seq, ("X",), tuple(events))


def pred(name="p1", seq="A", events=()):
    return PredictedAllele(name, seq, ("X",), 1, -8, tuple(events))


def evaluate(haps, preds, status="completed"):
    return evaluate_sample(TruthSample("s", tuple(haps)), RunObservation(status, tuple(preds)))


def test_swapped_sequence_assignment():
    result = match_alleles((hap("h1", "A"), hap("h2", "C")), (pred("p1", "C"), pred("p2", "A")))
    assert result.alternatives == (((0, 1), (1, 0)),)
    assert result.cost == 0


def test_canonical_count_consistency_uses_selected_reference_layout():
    prediction = PredictedAllele("p", "A", ("X",), 1, -1, fixed_repeat_count=2)
    result = evaluate([hap()], [prediction])
    pair = result["alternatives"][0]["pairs"][0]
    assert pair["fixed_repeat_count"] == 2
    assert pair["canonical_plus_fixed_matches_reported"] is True
    assert pair["canonical_plus9_matches_reported"] is False  # Legacy diagnostic only.


def test_wrong_site_parent_name_and_duplicate_calls():
    expected = [Event(25, "X", "dupC")]
    for wrong in [Event(999, "X", "dupC"), Event(25, "A", "dupC"), Event(25, "X", "dupCCCC")]:
        assert score_events(expected, [wrong]) == (0, 1, 1)
    assert score_events(expected, expected * 2) == (1, 0, 1)


def test_extra_allele_is_not_exact_and_events_are_false_positive():
    result = evaluate(
        [hap("h1", "A"), hap("h2", "C")],
        [pred("p1", "A"), pred("p2", "C"), pred("p3", "G", [Event(1, "X", "dupC")])],
    )
    assert result["metrics"]["all_sequences_exact"] == {"min": 0, "max": 0}
    assert result["metrics"]["event_fp"] == {"min": 1, "max": 1}
    assert result["extra_alleles"] == 1


def test_extra_allele_is_not_a_resolved_reconstruction_call():
    from dataclasses import replace

    predictions = [
        replace(pred("p1", "A"), sequence_source="p1"),
        replace(pred("p2", "C"), sequence_source="p2"),
        replace(pred("p3", "G"), sequence_source="p3"),
    ]
    result = evaluate([hap("h1", "A"), hap("h2", "C")], predictions)
    assert result["missing_alleles"] == 0
    assert result["extra_alleles"] == 1
    assert aggregate([result])["totals"]["resolved_reconstruction_call_rate"] == {
        "numerator": 0,
        "denominator": 1,
        "value": 0.0,
    }


def test_assignment_ties_use_conservative_event_bounds():
    event = Event(1, "X", "dupC")
    result = evaluate(
        [hap("h1", "A", [event]), hap("h2", "A")], [pred("p1", "A", [event]), pred("p2", "A")]
    )
    assert result["assignment_ties"] == 2
    assert result["metrics"]["event_tp"] == {"min": 0, "max": 1}
    assert result["metrics"]["event_fn"] == {"min": 0, "max": 1}
    assert result["metrics"]["sequence_exact"] == {"min": 2, "max": 2}
    assert result["metrics"]["all_sequences_exact"] == {"min": 0, "max": 0}


def test_failed_stale_predictions_keep_denominator_and_no_normal_tn():
    event = Event(1, "X", "dupC")
    failed = evaluate([hap(events=[event]), hap("h2")], [pred(events=[event])], "execution_failed")
    assert failed["metrics"]["event_fn"]["min"] == 1
    assert failed["missing_alleles"] == 2
    normal = evaluate([hap()], [], "not_attempted")
    assert normal["normal_true_negative"] is False
    totals = aggregate([failed, normal])["totals"]
    assert totals["truth_haplotypes"] == 3
    assert totals["event_recall"]["denominator"] == 1
    assert totals["event_recall"]["value"] == 0


def test_empty_and_unsupported_candidates():
    event = Event(1, "X", "dupC")
    result = evaluate([hap(events=[event])], [pred(events=[event])])
    assert result["metrics"]["event_tp"]["min"] == 1
    assert result["metrics"]["supported_event_tp"]["min"] == 0
    assert aggregate([])["totals"]["event_recall"]["value"] is None
    result = evaluate([hap()], [])
    assert result["metrics"]["all_sequences_exact"]["min"] == 0


def test_missing_wildtype_is_reconstruction_loss_not_mutation_fn():
    result = evaluate([hap("h1", "A"), hap("h2", "C")], [pred()])
    assert result["missing_alleles"] == 1
    assert result["metrics"]["event_fn"]["max"] == 0
    assert result["metrics"]["all_counts_exact"]["max"] == 0


def test_haplotype_two_mutations_cannot_match_haplotype_one():
    event = Event(1, "X", "dupC", True, True, True, "exact_sequence_concordance")
    result = evaluate(
        [hap("h1", "A"), hap("h2", "C", [event])], [pred("p1", "A", [event]), pred("p2", "C")]
    )
    assert result["metrics"]["event_tp"]["max"] == 0
    assert result["metrics"]["event_fn"]["min"] == 1
    assert result["metrics"]["event_fp"]["min"] == 1
    assert result["metrics"]["supported_event_fp"]["min"] == 1


def test_literal_ambiguity_and_count_error():
    result = evaluate([hap(seq="A")], [PredictedAllele("p", "N", ("X",), 3, -6)])
    detail = result["alternatives"][0]["pairs"][0]
    assert detail["sequence_edit_distance"] == 1
    assert detail["ambiguous_bases"] == 1
    assert detail["repeat_count_error"] == 2
    assert result["metrics"]["count_within1"]["min"] == 0
    assert result["metrics"]["count_within2"]["min"] == 1


def test_identical_consensus_copies_require_independent_genotype_evidence():
    from dataclasses import replace

    truth = [hap("h1", "A"), hap("h2", "A")]
    copies = [pred("p1", "A"), pred("p2", "A")]
    result = evaluate(truth, copies)
    assert result["metrics"]["all_sequences_exact"]["max"] == 0
    assert result["metrics"]["sequence_exact"] == {"min": 2, "max": 2}
    assert result["metrics"]["independent_sequence_exact"] == {"min": 1, "max": 1}
    assert result["missing_alleles"] == 0  # Literal output cardinality remains explicit.
    assert result["unproven_duplicate_alleles"] == 1
    assert result["independent_missing_alleles"] == 1
    assert any("Unproven duplicate" in warning for warning in result["warnings"])
    totals = aggregate([result])["totals"]
    assert totals["independent_missing_alleles"] == 1
    assert totals["literal_sequence_accuracy"]["value"] == 1.0
    assert totals["sequence_accuracy"]["value"] == 0.5
    assert totals["literal_structure_accuracy"]["value"] == 1.0
    assert totals["structure_accuracy"]["value"] == 0.5
    assert totals["independently_recovered_sequence_accuracy"] == totals["sequence_accuracy"]
    first = replace(copies[0], sequence_source="block:GT1", independent_haplotype_evidence=True)
    second = replace(copies[1], sequence_source="block:GT2", independent_haplotype_evidence=True)
    evidenced = evaluate(truth, [first, second])
    assert evidenced["metrics"]["all_sequences_exact"]["min"] == 1
    assert evidenced["metrics"]["all_counts_exact"]["min"] == 1
    assert evidenced["metrics"]["independent_sequence_exact"] == {"min": 2, "max": 2}
    assert evidenced["independent_missing_alleles"] == 0
    assert evidenced["unproven_duplicate_alleles"] == 0
    assert (
        evaluate(truth, [first, replace(second, sequence_source="block:GT1")])["metrics"][
            "all_sequences_exact"
        ]["max"]
        == 0
    )
    empty_source = replace(first, sequence_source="")
    empty_source_result = evaluate(truth, [empty_source, second])
    assert empty_source_result["metrics"]["independent_sequence_exact"] == {
        "min": 1,
        "max": 1,
    }


def test_unproven_duplicate_structure_credit_retains_assignment_bounds():
    first = pred("p1", "A")
    second = PredictedAllele("p2", "A", ("Y",), 1, -8)
    truth = [hap("h1", "A"), TruthHaplotype("h2", "A", ("Y",))]
    result = evaluate(truth, [first, second])
    assert result["assignment_ties"] == 2
    assert result["metrics"]["structure_exact"] == {"min": 0, "max": 2}
    assert result["metrics"]["independent_structure_exact"] == {"min": 0, "max": 1}


def test_sequence_exactness_does_not_depend_on_scorer(monkeypatch):
    monkeypatch.setattr("muc_one_span.evaluation.matching.edit_distance", lambda a, b: 0)
    result = evaluate([hap(seq="A")], [pred(seq="C")])
    assert result["metrics"]["sequence_exact"]["max"] == 0
    assert result["alternatives"][0]["pairs"][0]["sequence_exact"] is False


def test_reused_source_does_not_establish_two_haplotypes():
    from dataclasses import replace

    first = replace(pred("p1", "A"), sequence_source="same-consensus-file")
    second = replace(pred("p2", "C"), sequence_source="same-consensus-file")
    result = evaluate([hap("h1", "A"), hap("h2", "C")], [first, second])
    assert result["metrics"]["all_sequences_exact"]["max"] == 0
    assert result["metrics"]["sequence_exact"] == {"min": 2, "max": 2}
    assert result["metrics"]["independent_sequence_exact"] == {"min": 1, "max": 1}


def test_supported_sample_detection_and_no_call_denominators():
    normal = evaluate([hap()], [pred(events=[Event(1, "X", "novel")])])
    failed = evaluate([hap()], [], "execution_failed")
    totals = aggregate([normal, failed])["totals"]
    assert totals["supported_normal_specificity"] == {
        "numerator": 1,
        "denominator": 1,
        "value": 1.0,
    }
    assert totals["normal_specificity"]["numerator"] == 0
    assert totals["no_call_rate"] == {"numerator": 1, "denominator": 2, "value": 0.5}


def test_exact_supported_metrics_exclude_legacy_proximity_boolean():
    expected = Event(1, "X", "dupC")
    legacy = Event(1, "X", "dupC", True, True, True, "unknown")
    exact = Event(1, "X", "dupC", True, True, True, "exact_sequence_concordance")

    assert not legacy.supported
    assert legacy.legacy_supported
    assert exact.supported
    legacy_result = evaluate([hap(events=[expected])], [pred(events=[legacy])])
    assert legacy_result["metrics"]["supported_event_tp"] == {"min": 0, "max": 0}
    assert legacy_result["metrics"]["supported_event_fn"] == {"min": 1, "max": 1}
    assert legacy_result["metrics"]["legacy_supported_event_tp"] == {"min": 1, "max": 1}
    totals = aggregate([legacy_result])["totals"]
    assert totals["supported_sample_alarm_sensitivity"]["value"] == 0.0
    assert totals["legacy_supported_sample_alarm_sensitivity"]["value"] == 1.0


def test_unmatched_legacy_supported_event_remains_a_legacy_false_positive():
    legacy = Event(1, "X", "dupC", True, True, True, "unknown")
    result = evaluate([hap("h1", "A")], [pred("p1", "A"), pred("extra", "C", [legacy])])
    assert result["metrics"]["supported_event_fp"] == {"min": 0, "max": 0}
    assert result["metrics"]["legacy_supported_event_fp"] == {"min": 1, "max": 1}


def test_supported_policy_text_matches_event_supported_semantics():
    from muc_one_span.evaluation.scoring import SCORING_POLICY

    assert SCORING_POLICY["supported_events"] == (
        "frameshift AND template_match AND vcf_support are true AND "
        "vcf_support_status is exact_sequence_concordance; OR frameshift AND "
        "template_match AND read_support_status is supported (hybrid engine)"
    )
    cases = (
        (True, True, True, "exact_sequence_concordance", "unknown"),  # legacy-path support
        (True, True, False, "exact_sequence_concordance", "unknown"),  # no vcf_support boolean
        (True, True, True, "unknown", "unknown"),  # status not exact concordance
        (True, True, False, "not_applicable_read_consensus", "supported"),  # read-only path
        (True, False, False, "unknown", "supported"),  # missing template_match
        (False, True, False, "unknown", "supported"),  # missing frameshift
        (True, True, False, "unknown", "discordant"),  # read support not "supported"
    )
    for frameshift, template_match, vcf_support, support_status, read_support_status in cases:
        event = Event(
            1,
            "X",
            "dupC",
            frameshift,
            template_match,
            vcf_support,
            support_status,
            read_support_status,
        )
        expected = (
            frameshift
            and template_match
            and vcf_support
            and support_status == "exact_sequence_concordance"
        ) or (frameshift and template_match and read_support_status == "supported")
        assert event.supported == expected


def test_ambiguous_positive_is_not_removed_from_specificity_denominator():
    normal = evaluate([hap()], [pred()])
    ambiguous_positive = evaluate(
        [hap()], [pred(events=[Event(1, "X", "dupC")])], "ambiguous_reconstruction"
    )
    ambiguous_negative = evaluate([hap()], [pred()], "ambiguous_reconstruction")
    totals = aggregate([normal, ambiguous_positive, ambiguous_negative])["totals"]
    assert totals["normal_specificity"] == {"numerator": 1, "denominator": 2, "value": 0.5}
    assert totals["normal_confident_negative_rate"]["denominator"] == 3
    assert totals["normal_unresolved_negative_count"] == 1
