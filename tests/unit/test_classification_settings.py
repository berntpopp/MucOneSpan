"""Observable classification and confidence runtime-setting contracts."""

from dataclasses import replace

import pytest

from muc_one_span import classify
from muc_one_span.classification_summary import _qual_to_confidence
from muc_one_span.classify import classify_repeat, classify_sequence, validate_mutations_against_vcf
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.settings import DEFAULT_SETTINGS


def test_novel_repeat_distance_setting_changes_only_the_label() -> None:
    rd = load_repeat_dictionary()
    sequence = "T" + rd.repeats["X"][1:]
    default = classify_repeat(sequence, rd)
    configured = classify_repeat(
        sequence,
        rd,
        settings=replace(DEFAULT_SETTINGS.classification, novel_repeat_edit_distance=0),
    )
    assert default["classification"] == "variant"
    assert configured["classification"] == "novel_repeat"
    assert configured["edit_distance"] == default["edit_distance"] == 1


def test_fit_distance_and_strict_setting_control_unresolved_segmentation() -> None:
    rd = load_repeat_dictionary()
    sequence = "T" + rd.repeats["X"][1:]
    settings = replace(
        DEFAULT_SETTINGS.classification,
        max_fit_edit_distance=0,
        strict_segmentation=True,
    )
    configured = classify_sequence(sequence, rd, settings=settings)
    overridden = classify_sequence(sequence, rd, settings=settings, strict_segmentation=False)
    assert configured["classified_bases"] == 0
    assert configured["reconstruction_status"] == "ambiguous_reconstruction"
    assert overridden["classified_bases"] == len(sequence)
    assert overridden["segmentation_policy"] == "candidate_windows"


def test_minimum_unit_fraction_controls_short_terminal_processing() -> None:
    rd = load_repeat_dictionary()
    sequence = rd.repeats["X"][:30]
    default = classify_sequence(sequence, rd)
    configured = classify_sequence(
        sequence,
        rd,
        settings=replace(DEFAULT_SETTINGS.classification, minimum_unit_fraction=0.6),
    )
    assert default["classified_bases"] == 30
    assert configured["classified_bases"] == 0


def test_max_indel_probe_controls_exact_variable_length_template() -> None:
    rd = load_repeat_dictionary()
    mutant = next(sequence for sequence in rd.mutated_sequences if len(sequence) != 60)
    default = classify_sequence(mutant, rd)
    configured = classify_sequence(
        mutant,
        rd,
        settings=replace(DEFAULT_SETTINGS.classification, max_indel_probe=0),
    )
    assert default["mutations_detected"][0]["template_match"] is True
    assert configured["structure"] != default["structure"]


def test_early_stop_distance_controls_alternate_probe_search(monkeypatch) -> None:
    rd = load_repeat_dictionary()
    sequence = "A" * 61

    def distance_by_size(window, _repeat_dict, *, settings=None):
        distance = {60: 2, 30: 5, 31: 1}.get(len(window), 5)
        return {
            "type": "unknown",
            "match": "closest",
            "closest_match": "X",
            "edit_distance": distance,
            "confidence": 0.5,
            "differences": [],
            "classification": "variant",
        }

    monkeypatch.setattr(classify, "classify_repeat", distance_by_size)
    default = classify_sequence(sequence, rd)
    configured = classify_sequence(
        sequence,
        rd,
        settings=replace(DEFAULT_SETTINGS.classification, early_stop_edit_distance=2),
    )
    assert default["repeats"][0]["edit_distance"] == 1
    assert configured["repeats"][0]["edit_distance"] == 2


def test_confidence_breakpoints_and_weights_are_configurable() -> None:
    settings = replace(
        DEFAULT_SETTINGS.confidence,
        qual_low=10.0,
        qual_high=30.0,
        weight_below=0.2,
        weight_low=0.4,
        weight_high=0.8,
    )
    assert _qual_to_confidence(9.0, settings=settings) == 0.2
    assert _qual_to_confidence(10.0, settings=settings) == 0.4
    assert _qual_to_confidence(20.0, settings=settings) == pytest.approx(0.6)
    assert _qual_to_confidence(30.0, settings=settings) == 0.8


def test_absent_and_boundary_weights_are_configurable_and_explicit_override_wins(
    monkeypatch,
) -> None:
    classification = {
        "mutations_detected": [{"repeat_index": 4}],
        "repeats": [{"confidence": 1.0} for _ in range(4)],
    }
    monkeypatch.setattr(
        "muc_one_span.variant_support.mutation_concordance",
        lambda *_: {4: ("absent", [])},
    )
    settings = replace(
        DEFAULT_SETTINGS.confidence,
        absent_weight=0.2,
        boundary_repeats=1,
        boundary_penalty=0.25,
    )
    configured = validate_mutations_against_vcf(classification, [], settings=settings)
    count_overridden = validate_mutations_against_vcf(
        classification,
        [],
        boundary_repeats=0,
        settings=settings,
    )
    both_overridden = validate_mutations_against_vcf(
        classification,
        [],
        boundary_repeats=1,
        boundary_penalty=0.5,
        settings=settings,
    )
    assert configured["repeats"][3]["confidence"] == 0.05
    assert count_overridden["repeats"][3]["confidence"] == 0.2
    assert both_overridden["repeats"][3]["confidence"] == 0.1
