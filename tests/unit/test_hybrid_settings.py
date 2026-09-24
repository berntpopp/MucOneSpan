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
