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
    assert (h.rejected_peak_noise_reads, h.hp_min_strand_reads) == (2, 5)
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


def test_hybrid_poa_and_polish_tunables_default_unchanged() -> None:
    # Task 6 (S3/S7): these fields replace literals formerly hardcoded in hybrid/poa.py
    # and hybrid/polish.py (POA sampling window, insertion-vote majority fraction, and
    # homopolymer minimum run length); their defaults reproduce the old behaviour.
    h = DEFAULT_SETTINGS.hybrid
    assert (h.poa_sample_window_floor_bp, h.poa_sample_window_frac) == (15.0, 0.006)
    assert (h.polish_insertion_majority_frac, h.hp_vote_min_run) == (0.5, 4)


def test_hybrid_assign_flank_bp_default_unchanged() -> None:
    # Task 7 (S5/S6): replaces the brief's bare module-level FLANK_BP = 500 literal with
    # a validated HybridSettings field; the default reproduces the brief's value.
    assert DEFAULT_SETTINGS.hybrid.assign_flank_bp == 500


def test_hybrid_engine_orchestration_defaults() -> None:
    # Task 11: replace the brief's n_poa * 3 polishing cap, MAX_QC_READS = 200 and the
    # one-unit partial-fragment floor with validated fields; defaults keep those values.
    h = DEFAULT_SETTINGS.hybrid
    assert (h.polish_max_reads, h.qc_residual_max_reads) == (3 * h.n_poa, 200)
    assert h.polish_partial_min_units == 1.0


def test_hybrid_event_evidence_defaults() -> None:
    # Task 9 (S8/S10): replace the brief's run >= 4, MAXLEN = 16, 0.5 pseudo-count and
    # residual run >= 3 literals with validated fields; defaults keep the brief's values.
    h = DEFAULT_SETTINGS.hybrid
    assert (h.hp_event_min_run, h.hp_max_run_len) == (4, 16)
    assert (h.hp_background_pseudocount, h.qc_residual_min_run) == (0.5, 3)


def test_hybrid_smear_significance_defaults() -> None:
    # Fix round 4 (C4.2): the smear decision is a one-sided significance test of the
    # candidate's core count against its local below-top background; these defaults
    # come from the fix-round-4 simulation sweep (see task-5-report.md). The round-2/3
    # excess-fraction fields and smear_background_floor are removed.
    h = DEFAULT_SETTINGS.hybrid
    assert (h.smear_test_alpha, h.smear_test_borderline_factor) == (0.001, 3.0)
    assert (h.smear_test_correction, h.smear_test_window_frac) == ("bonferroni", 0.25)
    assert (h.smear_background_flank_units, h.smear_background_min_reads) == (2.0, 5)
    for removed in (
        "smear_min_expected",
        "smear_explained_frac",
        "smear_confident_frac",
        "smear_low_background_min_support",
        "smear_background_floor",
    ):
        assert not hasattr(h, removed)


def test_hybrid_phase_tunables_default_unchanged() -> None:
    # Task 8 (S4): replace the literals of the prototype linked-site split (hetsplit.py):
    # site-read cap 300, run min length 3, background window +-3, minor floor 5 reads,
    # run background multiplier 4 (D6), gap AF factor 1.5 and a 10-read floor for
    # pairwise linkage. Fix round 1 removed phase_max_minor_alleles: run sites scan all
    # alleles (as the prototype does) and other sites use only the top minor.
    h = DEFAULT_SETTINGS.hybrid
    assert (h.phase_max_site_reads, h.phase_run_min_len, h.phase_run_bg_window) == (300, 3, 3)
    assert h.phase_min_minor_reads == 5 and not hasattr(h, "phase_max_minor_alleles")
    assert (h.phase_run_bg_multiplier, h.phase_gap_af_factor) == (4.0, 1.5)
    assert h.phase_min_pair_reads == 10
    # Fix round 2: strand consistency is a strand-bias test at this alpha.
    assert h.phase_strand_bias_alpha == 0.001


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("phase_max_site_reads", 0),
        ("phase_run_min_len", 1),
        ("phase_run_bg_window", 0),
        ("phase_min_minor_reads", 0),
        ("phase_run_bg_multiplier", -0.5),
        ("phase_gap_af_factor", 0.5),
        ("phase_min_pair_reads", 1),
        ("phase_strand_bias_alpha", 0.0),
        ("phase_strand_bias_alpha", 1.0),
        ("phase_max_site_reads", 1.5),
        ("n_poa", 0),
        ("assign_margin", -1),
        ("het_af_min", 0.9),
        ("poa_backend", "medaka"),
        ("assign_max_error_rate", 1.5),
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
        ("smear_test_alpha", 0.0),
        ("smear_test_alpha", 1.0),
        ("smear_test_borderline_factor", 0.5),
        ("smear_test_correction", "holm"),
        ("smear_test_window_frac", 0.0),
        ("smear_test_window_frac", 1.5),
        ("smear_background_flank_units", 0.0),
        ("smear_background_min_reads", 0),
        ("poa_sample_window_floor_bp", -1.0),
        ("poa_sample_window_frac", -0.1),
        ("polish_insertion_majority_frac", -0.1),
        ("polish_insertion_majority_frac", 1.5),
        ("hp_vote_min_run", 1),
        ("hp_vote_min_run", -1),
        ("assign_flank_bp", 0),
        ("assign_flank_bp", -1),
        ("hp_event_min_run", 1),
        ("hp_max_run_len", 4),
        ("hp_background_pseudocount", 0.0),
        ("hp_background_pseudocount", -0.5),
        ("qc_residual_min_run", 1),
        ("hp_min_reads", 0),
        ("hp_llr_min", 0.0),
        ("hp_llr_min", -1.0),
        ("polish_max_reads", 0),
        ("qc_residual_max_reads", 0),
        ("polish_partial_min_units", -1.0),
    ],
)
def test_hybrid_rejects_invalid(field: str, value: object) -> None:
    with pytest.raises(ValueError, match=f"hybrid.{field}"):
        HybridSettings(**{field: value})  # type: ignore[arg-type]


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
