"""Hybrid-engine settings: defaults, validation, and configuration round trip."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from muc_one_span.settings import (
    DEFAULT_SETTINGS,
    HybridSettings,
    RunSettings,
    load_settings,
    settings_as_dict,
)


def test_hybrid_defaults() -> None:
    h = DEFAULT_SETTINGS.hybrid
    assert (h.anchor_max_edits, h.n_poa, h.assign_margin, h.depth_adequate_spanning) == (
        12,
        40,
        3,
        30,
    )
    assert (h.poa_backend, h.assign_max_error_rate, h.max_unassigned_spanning_fraction) == (
        "pyabpoa",
        0.15,
        0.2,
    )
    assert (h.smear_min_prominence, h.rejected_peak_noise_reads, h.hp_min_strand_reads) == (
        3.0,
        2,
        5,
    )
    assert (DEFAULT_SETTINGS.run.engine, DEFAULT_SETTINGS.run.assay) == ("ladder", "amplicon")


def test_hybrid_length_model_and_anchor_tunables_default_unchanged() -> None:
    # Fix round 1: these fields replace literals formerly hardcoded in hybrid/spans.py and
    # hybrid/lengths.py; their defaults reproduce the old, already-validated behaviour.
    h = DEFAULT_SETTINGS.hybrid
    assert (h.flank_anchor_bp, h.flank_anchor_edit_divisor, h.flank_anchor_edit_floor) == (
        30,
        4,
        2,
    )
    assert (h.kde_bandwidth_base_bp, h.kde_bandwidth_per_bp, h.kde_kernel_truncation_bw) == (
        8.0,
        0.004,
        4.0,
    )
    assert (h.kde_grid_step_bp, h.kde_grid_margin_bp) == (2.0, 100.0)
    assert (h.smear_short_product_units, h.peak_far_near_boundary_units) == (1.5, 2.0)
    assert h.peak_min_separation_units == 0.7
    assert h.smear_background_floor == 1.0


def test_hybrid_smear_model_defaults_unchanged() -> None:
    # Fix round 2 (C4.2): one explicit smear model replaces the round-1 shoulder/ratio
    # fields (smear_shoulder_width_mult, smear_shoulder_floor, smear_background_ratio_min,
    # now removed) with these four; defaults are the tuned values from the fix-round-2
    # simulation sweep (see task-5-report.md).
    h = DEFAULT_SETTINGS.hybrid
    assert (h.smear_min_expected, h.smear_explained_frac, h.smear_confident_frac) == (
        0.9,
        0.10,
        0.13,
    )
    assert h.smear_low_background_min_support == 5


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("n_poa", 0),
        ("assign_margin", -1),
        ("het_af_min", 0.9),
        ("poa_backend", "medaka"),
        ("assign_max_error_rate", 1.5),
        ("smear_min_prominence", 0.5),
        ("max_span_units", 10),
        ("depth_adequate_spanning", 5),
        ("hp_vote", 1),
        ("flank_anchor_bp", 0),
        ("flank_anchor_edit_divisor", 0),
        ("flank_anchor_edit_floor", -1),
        ("kde_bandwidth_base_bp", 0.0),
        ("kde_bandwidth_per_bp", -0.1),
        ("kde_kernel_truncation_bw", 0.5),
        ("kde_grid_step_bp", 0.0),
        ("kde_grid_margin_bp", -1.0),
        ("smear_short_product_units", 0.0),
        ("peak_far_near_boundary_units", -1.0),
        ("peak_min_separation_units", -0.1),
        ("smear_background_floor", 0.0),
        ("smear_min_expected", -1.0),
        ("smear_explained_frac", -0.1),
        ("smear_confident_frac", 1.5),
        ("smear_low_background_min_support", -1),
    ],
)
def test_hybrid_rejects_invalid(field: str, value: object) -> None:
    with pytest.raises(ValueError, match=f"hybrid.{field}"):
        HybridSettings(**{field: value})  # type: ignore[arg-type]


def test_hybrid_rejects_confident_frac_below_explained_frac() -> None:
    with pytest.raises(ValueError, match=r"hybrid\.smear_confident_frac"):
        HybridSettings(smear_explained_frac=0.2, smear_confident_frac=0.1)


@pytest.mark.parametrize(("field", "value"), [("engine", "assembly"), ("assay", "wgs")])
def test_run_engine_and_assay_are_validated(field: str, value: str) -> None:
    with pytest.raises(ValueError, match=f"run.{field}"):
        RunSettings(**{field: value})  # type: ignore[arg-type]


def test_hybrid_roundtrip(tmp_path: Path) -> None:
    data = settings_as_dict(DEFAULT_SETTINGS)
    data["hybrid"]["n_poa"] = 25
    data["run"]["engine"] = "hybrid"
    path = tmp_path / "s.json"
    path.write_text(json.dumps(data))
    loaded = load_settings(path)
    assert loaded.hybrid.n_poa == 25 and loaded.run.engine == "hybrid"


def test_partial_hybrid_section_keeps_other_defaults(tmp_path: Path) -> None:
    path = tmp_path / "s.json"
    path.write_text('{"schema_version": 1, "hybrid": {"poa_backend": "pyspoa"}}')
    loaded = load_settings(path)
    assert loaded.hybrid.poa_backend == "pyspoa" and loaded.hybrid.n_poa == 40
