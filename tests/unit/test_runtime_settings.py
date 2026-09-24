"""Strict runtime settings preserve defaults and reject invalid scientific inputs."""

import json
import logging
import warnings
from dataclasses import FrozenInstanceError
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.settings import (
    DEFAULT_SETTINGS,
    AlleleSelectionSettings,
    CallingSettings,
    ClassificationSettings,
    ConfidenceSettings,
    ConsensusSettings,
    ReadPhasingSettings,
    ReferenceLayoutSettings,
    RunSettings,
    RuntimeSettings,
    load_settings,
    settings_as_dict,
)


def test_default_configuration_roundtrip_and_immutability(tmp_path: Path) -> None:
    settings = load_settings(None)
    assert settings == DEFAULT_SETTINGS == RuntimeSettings()
    assert settings.run.threads == 4 and settings.run.min_coverage == 10
    assert settings.run.min_qual == 5 and settings.run.platform == "hifi"
    assert settings.classification.max_indel_probe == 30
    assert settings.classification.minimum_unit_fraction == 0.5
    assert settings.confidence.qual_low == 5 and settings.confidence.qual_high == 20
    assert settings.read_phasing.internal_downsampling is None
    assert settings.calling.read_phase is False
    assert settings.reference_layout.fixed_repeat_count == 9
    assert settings.reference_layout.left_anchor_id == "1"
    assert settings.reference_layout.right_anchor_id == "9"
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(settings_as_dict(settings)))
    assert load_settings(path) == settings
    with pytest.raises(FrozenInstanceError):
        settings.run.threads = 8  # type: ignore[misc]


@pytest.mark.parametrize(
    "contents",
    [
        "[]",
        "null",
        "{}",
        '{"schema_version":2}',
        '{"schema_version":true}',
        '{"schema_version":1.0}',
        '{"schema_version":1,"unknown":{}}',
        '{"schema_version":1,"run":null}',
        '{"schema_version":1,"run":[]}',
        '{"schema_version":1,"run":{"unknown":2}}',
        '{"schema_version":1,"schema_version":1}',
        '{"schema_version":1,"run":{"threads":2,"threads":4}}',
        '{"schema_version":1,"run":{"min_qual":NaN}}',
        '{"schema_version":1,"run":{"min_qual":Infinity}}',
        '{"schema_version":1,"run":{"min_qual":1e999}}',
        '{"schema_version":1,"repeat_dictionary":false}',
        '{"schema_version":1,"run":{"threads":',
    ],
)
def test_bad_json_configuration_is_rejected(tmp_path: Path, contents: str) -> None:
    path = tmp_path / "bad.json"
    path.write_text(contents)
    with pytest.raises(ValueError):
        load_settings(path)


@pytest.mark.parametrize(
    ("section", "values"),
    [
        (RunSettings, {"threads": 0}),
        (RunSettings, {"threads": True}),
        (RunSettings, {"threads": 1.0}),
        (RunSettings, {"min_coverage": -1}),
        (RunSettings, {"min_coverage": 0}),
        (RunSettings, {"min_qual": -1}),
        (RunSettings, {"min_qual": float("nan")}),
        (RunSettings, {"min_qual": "5"}),
        (RunSettings, {"platform": "pacbio"}),
        (RunSettings, {"report": 1}),
        (RunSettings, {"reference": ""}),
        (RunSettings, {"clair3_model": None}),
        (RunSettings, {"minimap2_preset": " "}),
        (AlleleSelectionSettings, {"min_gap": -1}),
        (AlleleSelectionSettings, {"min_gap": 0}),
        (AlleleSelectionSettings, {"valley_min_points": 2}),
        (AlleleSelectionSettings, {"valley_min_separation": 0}),
        (AlleleSelectionSettings, {"refinement_max_shift": -1}),
        (ClassificationSettings, {"max_indel_probe": -1}),
        (ClassificationSettings, {"max_fit_edit_distance": -1}),
        (ClassificationSettings, {"novel_repeat_edit_distance": -1}),
        (ClassificationSettings, {"minimum_unit_fraction": 0}),
        (ClassificationSettings, {"minimum_unit_fraction": 1.1}),
        (ClassificationSettings, {"early_stop_edit_distance": -1}),
        (ClassificationSettings, {"strict_segmentation": "false"}),
        (ConsensusSettings, {"flank_length": -1}),
        (ConsensusSettings, {"anchor_bases": 0}),
        (ConsensusSettings, {"anchor_tolerance": -1}),
        (ConfidenceSettings, {"qual_low": -1}),
        (ConfidenceSettings, {"qual_low": 20}),
        (ConfidenceSettings, {"qual_high": float("inf")}),
        (ConfidenceSettings, {"weight_below": -0.1}),
        (ConfidenceSettings, {"weight_low": 1.1}),
        (ConfidenceSettings, {"weight_high": True}),
        (ConfidenceSettings, {"absent_weight": -0.1}),
        (ConfidenceSettings, {"boundary_repeats": -1}),
        (ConfidenceSettings, {"boundary_penalty": 2}),
        (CallingSettings, {"sample_name": ""}),
        (CallingSettings, {"sample_name": "two samples"}),
        (CallingSettings, {"sample_name": "a\tb"}),
        (CallingSettings, {"sample_name": "a\nb"}),
        (CallingSettings, {"sample_name": "a\x00b"}),
        (CallingSettings, {"sample_name": "a\x7fb"}),
        (CallingSettings, {"read_phase": 1}),
        (ReadPhasingSettings, {"internal_downsampling": 0}),
        (ReadPhasingSettings, {"mapping_quality": -1}),
        (ReadPhasingSettings, {"mapping_quality": False}),
        (ReferenceLayoutSettings, {"pre": ()}),
        (ReferenceLayoutSettings, {"pre": ["1"]}),
        (ReferenceLayoutSettings, {"pre": ("1", "1")}),
        (ReferenceLayoutSettings, {"pre": ("6",)}),
        (ReferenceLayoutSettings, {"after": ("",)}),
        (ReferenceLayoutSettings, {"after": (9,)}),
        (RuntimeSettings, {"schema_version": 2}),
        (RuntimeSettings, {"run": {}}),
        (RuntimeSettings, {"repeat_dictionary": 1}),
    ],
)
def test_direct_construction_validates_types_and_ranges(
    section: Any, values: dict[str, Any]
) -> None:
    with pytest.raises(ValueError):
        section(**values)


def test_partial_sections_resolve_paths_and_keep_unspecified_defaults(tmp_path: Path) -> None:
    path = tmp_path / "config.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "run": {"reference": "refs/a.fa", "clair3_model": "models/hifi", "threads": 8},
                "repeat_dictionary": "dictionary.json",
                "reference_layout": {"pre": ["first", "second"], "after": ["last"]},
                "read_phasing": {"internal_downsampling": 20, "mapping_quality": 0},
            }
        )
    )
    settings = load_settings(path)
    assert settings.run.reference == str(tmp_path / "refs/a.fa")
    assert settings.run.clair3_model == str(tmp_path / "models/hifi")
    assert settings.repeat_dictionary == str(tmp_path / "dictionary.json")
    assert settings.run.threads == 8 and settings.run.min_coverage == 10
    assert settings.reference_layout.pre == ("first", "second")
    assert settings.reference_layout.fixed_repeat_count == 3
    assert settings.reference_layout.left_anchor_id == "first"
    assert settings.reference_layout.right_anchor_id == "last"
    settings.reference_layout.validate_repeats({"first": "AC", "second": "GT", "last": "AC"})
    with pytest.raises(ValueError, match="second"):
        settings.reference_layout.validate_repeats({"first": "AC", "last": "AC"})
    assert settings.read_phasing.mapping_quality == 0


def test_zero_thresholds_and_closed_weight_boundaries_are_valid() -> None:
    settings = RuntimeSettings(
        run=RunSettings(min_coverage=1, min_qual=0),
        classification=ClassificationSettings(
            max_indel_probe=0,
            max_fit_edit_distance=0,
            novel_repeat_edit_distance=0,
            early_stop_edit_distance=0,
            minimum_unit_fraction=1,
        ),
        confidence=ConfidenceSettings(
            qual_low=0,
            qual_high=1,
            weight_below=0,
            weight_low=0,
            weight_high=1,
            absent_weight=0,
            boundary_repeats=0,
            boundary_penalty=1,
        ),
        consensus=ConsensusSettings(flank_length=0, anchor_tolerance=0),
    )
    assert settings.classification.max_fit_edit_distance == 0
    assert settings.consensus.flank_length == 0


def test_nondefault_sections_roundtrip_without_losing_values(tmp_path: Path) -> None:
    settings = RuntimeSettings(
        run=RunSettings(
            threads=2,
            platform="ont",
            min_coverage=7,
            min_qual=9,
            clair3_model=str(tmp_path / "model"),
            reference=str(tmp_path / "ref.fa"),
            minimap2_preset="map-ont",
            report=True,
        ),
        allele_selection=AlleleSelectionSettings(2, 5, 4, 0),
        classification=ClassificationSettings(15, 6, 4, 0.75, 0, True),
        consensus=ConsensusSettings(200, 15, 25),
        confidence=ConfidenceSettings(4, 30, 0.1, 0.6, 0.9, 0.2, 2, 0.4),
        calling=CallingSettings("alternate_sample", True),
        read_phasing=ReadPhasingSettings(20, 10),
        reference_layout=ReferenceLayoutSettings(("first",), ("last",)),
        repeat_dictionary=str(tmp_path / "dictionary.json"),
    )
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(settings_as_dict(settings)))
    assert load_settings(path) == settings


@pytest.mark.parametrize("value", ['"1"', "1", "null", '["1", "1"]'])
def test_layout_json_requires_distinct_string_arrays(tmp_path: Path, value: str) -> None:
    path = tmp_path / "settings.json"
    path.write_text('{"schema_version":1,"reference_layout":{"pre":' + value + "}}")
    with pytest.raises(ValueError):
        load_settings(path)


def test_unrepresentable_numeric_threshold_is_a_configuration_error() -> None:
    with pytest.raises(ValueError, match="finite number"):
        RunSettings(min_qual=10**1000)


@pytest.mark.parametrize(
    ("section", "values"),
    [
        ("run", {"min_coverage": 0}),
        ("allele_selection", {"min_gap": 0}),
        ("calling", {"sample_name": "invalid sample"}),
        ("calling", {"sample_name": "invalid\tcolumn"}),
        ("calling", {"sample_name": "invalid\x01control"}),
        ("reference_layout", {"min_units": 0}),
        ("reference_layout", {"min_units": True}),
        ("reference_layout", {"min_units": 20, "max_units": 10}),
    ],
)
def test_unsafe_stage_values_rejected_from_json(
    tmp_path: Path, section: str, values: dict[str, Any]
) -> None:
    path = tmp_path / "settings.json"
    path.write_text(json.dumps({"schema_version": 1, section: values}))
    with pytest.raises(ValueError):
        load_settings(path)


def test_reference_layout_range_roundtrip_and_validation(tmp_path: Path) -> None:
    assert DEFAULT_SETTINGS.reference_layout.min_units == 1
    assert DEFAULT_SETTINGS.reference_layout.max_units == 150
    layout = ReferenceLayoutSettings(min_units=20, max_units=30)
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(settings_as_dict(RuntimeSettings(reference_layout=layout))))
    assert load_settings(path).reference_layout == layout
    assert layout.fixed_repeat_count == 9
    for values in (
        {"min_units": 0},
        {"min_units": True},
        {"max_units": 0},
        {"max_units": 1.0},
        {"min_units": 20, "max_units": 10},
    ):
        with pytest.raises(ValueError):
            ReferenceLayoutSettings(**values)


def test_haploid_fraction_settings_are_validated() -> None:
    calling = CallingSettings()
    assert (calling.haploid_alt_fraction, calling.haploid_ref_fraction) == (0.5, 0.2)
    with pytest.raises(ValueError, match="haploid_ref_fraction must be below"):
        CallingSettings(haploid_alt_fraction=0.3, haploid_ref_fraction=0.3)
    with pytest.raises(ValueError, match="finite number"):
        CallingSettings(haploid_alt_fraction=1.5)


def test_selection_gate_settings_are_validated() -> None:
    selection = AlleleSelectionSettings()
    assert (selection.secondary_mode_min_fraction, selection.min_allele_primary_records) == (
        0.2,
        30,
    )
    with pytest.raises(ValueError, match="secondary_mode_min_fraction must be > 0"):
        AlleleSelectionSettings(secondary_mode_min_fraction=0)
    with pytest.raises(ValueError, match="min_allele_primary_records"):
        AlleleSelectionSettings(min_allele_primary_records=0)


def test_null_haploid_min_qual_roundtrips(tmp_path: Path) -> None:
    path = tmp_path / "settings.json"
    path.write_text('{"schema_version": 1, "calling": {"haploid_min_qual": null}}')
    settings = load_settings(path)
    assert settings.calling.haploid_min_qual is None
    path.write_text(json.dumps(settings_as_dict(settings)))
    assert load_settings(path) == settings


@pytest.mark.parametrize("consensus", ['{"haploid_min_qual": 3.0}', '{"haploid_majority": false}'])
def test_deprecated_consensus_haploid_settings_are_logged_when_changed(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, consensus: str
) -> None:
    """The deprecation reaches CLI users: logged, not a filtered DeprecationWarning."""
    path = tmp_path / "settings.json"
    path.write_text(f'{{"schema_version": 1, "consensus": {consensus}}}')
    with caplog.at_level(logging.WARNING, logger="muc_one_span.settings"):
        load_settings(path)
    assert [r.levelno for r in caplog.records] == [logging.WARNING]
    assert "consensus.haploid" in caplog.text


def test_deprecated_consensus_haploid_settings_silent_at_defaults(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path = tmp_path / "settings.json"
    for config in (
        settings_as_dict(DEFAULT_SETTINGS),
        {"schema_version": 1, "consensus": {"flank_length": 200}},
    ):
        path.write_text(json.dumps(config))
        with caplog.at_level(logging.WARNING), warnings.catch_warnings():
            warnings.simplefilter("error")
            load_settings(path)
    assert "consensus.haploid" not in caplog.text
