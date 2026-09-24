"""Benchmark sets: named technical factor mixes (standard headline, clean control, stress)."""

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import (
    DEFAULT_BENCH_CONFIG,
    GENERATION_SECTIONS,
    load_bench_config,
)
from muc_one_span.benchsim.bench_sets import BenchSet, FactorLevels, SetsConfig
from muc_one_span.benchsim.design import PROFILES
from muc_one_span.benchsim.realism_targets import load_targets, target_section

SETS = DEFAULT_BENCH_CONFIG.sets
AMPLICON = ("ont_amplicon_r10", "hifi_amplicon")


def _write(tmp_path: Path, data: object) -> Path:
    path = tmp_path / "bench.json"
    path.write_text(json.dumps(data))
    return path


def _levels(**overrides: Any) -> dict[str, Any]:
    base = SETS.definitions["clean"].profiles["hifi_amplicon"]
    data = {k: list(v) for k, v in base.__dict__.items()}
    return data | overrides


def _set(**overrides: Any) -> dict[str, Any]:
    data: dict[str, Any] = {
        "description": "x",
        "offpeak_share_cap": None,
        "profiles": {p: _levels() for p in PROFILES},
    }
    return data | overrides


def test_standard_is_default_headline_and_listed_first() -> None:
    assert SETS.default == SETS.headline == "standard"
    assert list(SETS.definitions) == ["standard", "clean", "stress"]
    assert SETS.order() == ("standard", "clean", "stress")
    assert all(set(s.profiles) == set(PROFILES) for s in SETS.definitions.values())


def test_standard_levels_are_realistic() -> None:
    standard = SETS.definitions["standard"].profiles
    assert standard["ont_amplicon_r10"].depths == (500, 1000, 2000)
    assert standard["hifi_amplicon"].depths == (200, 500, 1000)
    assert standard["ont_genomic_targeted"].depths == (30, 60, 100)
    for levels in standard.values():
        assert levels.error_levels == ("calibrated",)
        assert set(levels.pcr_levels) <= {"none", "calibrated"}
    for profile in AMPLICON:
        assert set(standard[profile].pcr_levels) == {"none", "calibrated"}


def test_standard_offpeak_is_the_typical_real_level() -> None:
    real = target_section(load_targets(), "ont_amplicon_r10")["span_off_gt1unit_frac"]
    standard = SETS.definitions["standard"]
    assert standard.offpeak_share_cap == real["median"]
    for profile in AMPLICON:
        levels = standard.profiles[profile]
        offpeak = max(levels.smear_levels) + max(levels.chimera_levels)
        assert offpeak <= real["median"] < real["max"]


def test_clean_is_the_standard_top_depth_without_artefacts() -> None:
    clean, standard = SETS.definitions["clean"], SETS.definitions["standard"]
    assert clean.offpeak_share_cap == 0
    for profile in PROFILES:
        levels = clean.profiles[profile]
        assert levels.depths == (max(standard.profiles[profile].depths),)
        assert levels.pcr_levels == ("none",) and levels.error_levels == ("calibrated",)
        assert levels.smear_levels == (0,) and levels.chimera_levels == (0,)


def test_stress_keeps_the_harsh_mix() -> None:
    stress = SETS.definitions["stress"]
    assert stress.offpeak_share_cap is None
    real_max = target_section(load_targets(), "ont_amplicon_r10")["span_off_gt1unit_frac"]["max"]
    for profile in PROFILES:
        levels = stress.profiles[profile]
        assert min(levels.depths) < DEFAULT_BENCH_CONFIG.atlas.min_resolvable_depth
        assert "poor" in levels.error_levels and "strong" in levels.pcr_levels
    amp = stress.profiles["ont_amplicon_r10"]
    assert max(amp.smear_levels) + max(amp.chimera_levels) > real_max


def test_sets_shape_generation() -> None:
    assert "sets" in GENERATION_SECTIONS
    stress = SETS.definitions["stress"]
    changed = replace(
        DEFAULT_BENCH_CONFIG,
        sets=replace(
            SETS, definitions={**SETS.definitions, "stress": replace(stress, description="y")}
        ),
    )
    assert changed.generation_sha256() != DEFAULT_BENCH_CONFIG.generation_sha256()


def test_json_replaces_definitions(tmp_path: Path) -> None:
    data = {
        "schema_version": 1,
        "sets": {
            "definitions": {"only": _set()},
            "default": "only",
            "headline": "only",
            "legacy": "only",
        },
    }
    cfg = load_bench_config(_write(tmp_path, data))
    assert list(cfg.sets.definitions) == ["only"]
    assert isinstance(cfg.sets.definitions["only"], BenchSet)
    levels = cfg.sets.definitions["only"].profiles["hifi_amplicon"]
    assert isinstance(levels, FactorLevels) and isinstance(levels.depths, tuple)
    assert cfg.sets.order() == ("only",)


@pytest.mark.parametrize(
    ("sets", "match"),
    [
        ({"default": "nope"}, r"sets\.default"),
        ({"headline": "nope"}, r"sets\.headline"),
        ({"legacy": "nope"}, r"sets\.legacy"),
        ({"definitions": {}}, "at least one set"),
        ({"definitions": {"a-b": _set()}}, "set name"),
        ({"definitions": {"standard": {"description": "x"}}}, "missing"),
        ({"definitions": {"standard": _set(extra=1)}}, "unknown"),
        ({"definitions": {"standard": _set(profiles={})}}, "at least one profile"),
        ({"definitions": {"standard": _set(profiles={"hifi_amplicon": []})}}, "JSON object"),
        (
            {"definitions": {"standard": _set(profiles={"x": _levels(depths=[0])})}},
            r"standard\.profiles\.x\.depths",
        ),
        (
            {"definitions": {"standard": _set(profiles={"x": _levels(pcr_levels=["hot"])})}},
            "pcr_levels",
        ),
        (
            {"definitions": {"standard": _set(profiles={"x": _levels(error_levels=[])})}},
            "error_levels",
        ),
        (
            {"definitions": {"standard": _set(profiles={"x": _levels(smear_levels=[2])})}},
            "smear_levels",
        ),
        (
            {
                "definitions": {
                    "standard": _set(
                        offpeak_share_cap=0.1,
                        profiles={"x": _levels(smear_levels=[0.1], chimera_levels=[0.05])},
                    )
                }
            },
            "exceeds offpeak_share_cap",
        ),
        ({"definitions": {"standard": _set(offpeak_share_cap=2)}}, "offpeak_share_cap"),
        ({"definitions": {"standard": _set(description=3)}}, "description"),
    ],
)
def test_invalid_sets_are_rejected(tmp_path: Path, sets: dict[str, Any], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        load_bench_config(_write(tmp_path, {"schema_version": 1, "sets": sets}))


def test_order_puts_the_headline_first() -> None:
    cfg = replace(SETS, headline="clean")
    assert cfg.order() == ("clean", "standard", "stress")
    assert cfg.order(["stress", "zzz", "standard"]) == ("standard", "stress", "zzz")


def test_set_config_rejects_names_outside_definitions_directly() -> None:
    with pytest.raises(ValueError, match=r"sets\.headline"):
        SetsConfig(headline="nope")
