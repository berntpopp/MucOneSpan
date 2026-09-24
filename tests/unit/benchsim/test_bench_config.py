import json
from dataclasses import replace
from pathlib import Path

import pytest

from muc_one_span.benchsim.bench_config import (
    DEFAULT_BENCH_CONFIG,
    BenchConfig,
    DesignConfig,
    load_bench_config,
)
from muc_one_span.benchsim.realism_targets import load_targets, target_section

CFG = DEFAULT_BENCH_CONFIG


def _write(tmp_path: Path, data: object) -> Path:
    path = tmp_path / "bench.json"
    path.write_text(json.dumps(data))
    return path


def test_defaults_are_valid_and_hash_is_stable() -> None:
    assert load_bench_config(None) is CFG
    assert BenchConfig().sha256() == CFG.sha256()
    cap = CFG.design.offpeak_share_cap
    assert json.loads(json.dumps(CFG.to_dict()))["design"]["offpeak_share_cap"] == cap


def test_json_overlays_defaults_and_coerces_lists(tmp_path: Path) -> None:
    floor = CFG.amount.min_minor_share / 2
    path = _write(
        tmp_path,
        {
            "schema_version": 1,
            "amount": {"min_minor_share": floor},
            "design": {"depths": {"hifi_amplicon": [7], "ont_amplicon_r10": [7], "x": [1]}},
        },
    )
    cfg = load_bench_config(path)
    assert cfg.amount.min_minor_share == floor
    assert cfg.design.depths["hifi_amplicon"] == (7,)
    assert cfg.realism == CFG.realism and cfg.sha256() != CFG.sha256()


@pytest.mark.parametrize(
    ("data", "match"),
    [
        ({"schema_version": 1, "nope": {}}, "unknown bench config sections"),
        ({"schema_version": 1, "report": {"nope": 1}}, "unknown report fields"),
        ({"report": {}}, "schema_version"),
        ({"schema_version": 1, "report": {"alpha": 0}}, "report.alpha"),
        ({"schema_version": 1, "design": {"chimera_levels": 0.1}}, "JSON array"),
        ({"schema_version": 1, "structures": []}, "JSON object"),
        ({"schema_version": 1, "realism": {"offset_bin_min": 5}}, "offset_bin_min"),
        ({"schema_version": 1, "design": {"compositions": {"markov": 0.5}}}, "sum to 1"),
        ({"schema_version": 1, "design": {"pcr_levels": ["hot"]}}, "pcr_levels"),
        ({"schema_version": 1, "amount": {"pcr_slope_per_unit": {"none": 0}}}, "pcr_slope"),
        ({"schema_version": 1, "design": {"delta_ranges": {">20": [21, 500]}}}, "length range"),
        ({"schema_version": 1, "design": {"depths": {"hifi_amplicon": [0]}}}, "positive"),
        ({"schema_version": 1, "design": {"split_sizes": {}}}, "at least one split"),
        ({"schema_version": 1, "design": {"split_sizes": {"dev": 3}}}, "every split"),
        ([], "JSON object"),
        (
            {"schema_version": 1, "design": {"compositions": {"markov": 1.5, "rare_units": -0.5}}},
            r"design\.compositions\.markov",
        ),
        (
            {
                "schema_version": 1,
                "design": {"compositions": {"markov": 1.0, "rare_units": -0.0001}},
            },
            r"design\.compositions\.rare_units",
        ),
        ({"schema_version": 1, "design": {"compositions": {"bogus": 1.0}}}, "bogus"),
        ({"schema_version": 1, "design": {"compositions": {"markov": 0.0}}}, "positive"),
        (
            {"schema_version": 1, "design": {"compositions": {"markov": "x"}}},
            r"compositions\.markov",
        ),
        ({"schema_version": 1, "design": {"delta_ranges": {"weird": [1, 2]}}}, "weird"),
        ({"schema_version": 1, "design": {"delta_ranges": {"1": [1]}}}, r"delta_ranges\.1"),
        ({"schema_version": 1, "design": {"delta_ranges": {"0_identical": [1, 1]}}}, "0_identical"),
        ({"schema_version": 1, "design": {"chimera_levels": [[1]]}}, "chimera_levels"),
        ({"schema_version": 1, "design": {"pcr_levels": [["none"]]}}, "pcr_levels"),
        ({"schema_version": 1, "design": {"depths": {"hifi_amplicon": [[5]]}}}, "depths"),
        ({"schema_version": 1, "design": {"split_sizes": {"dev": "x"}}}, "split_sizes"),
        ({"schema_version": 1, "amount": {"pcr_slope_per_unit": {"none": "x"}}}, "pcr_slope"),
    ],
)
def test_invalid_config_is_rejected(tmp_path: Path, data: object, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        load_bench_config(_write(tmp_path, data))


def test_duplicate_keys_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "bench.json"
    path.write_text('{"schema_version": 1, "run": {"threads": 1, "threads": 2}}')
    with pytest.raises(ValueError, match="duplicate"):
        load_bench_config(path)


def test_regular_smear_levels_are_capped_stress_is_exempt() -> None:
    design = CFG.design
    worst_chimera = max(design.chimera_levels)
    for split, levels in design.smear_levels.items():
        offpeak = max(levels) + worst_chimera
        assert (offpeak <= design.offpeak_share_cap) == (split not in design.offpeak_cap_exempt)
    too_high = {**design.smear_levels, "dev": design.smear_levels["stress"]}
    with pytest.raises(ValueError, match="exceeds offpeak_share_cap"):
        replace(design, smear_levels=too_high)
    DesignConfig(smear_levels=too_high, offpeak_cap_exempt=("dev", "stress"))


def test_offpeak_cap_is_the_public_real_maximum() -> None:
    amplicon = target_section(load_targets(), "ont_amplicon_r10")
    assert CFG.design.offpeak_share_cap == amplicon["span_off_gt1unit_frac"]["max"]


def test_realism_bins_match_the_public_target_file() -> None:
    amplicon = target_section(load_targets(), "ont_amplicon_r10")
    bins = amplicon["span_offset_pmf_15bp_bins"]
    assert set(bins) == set(CFG.realism.size_keys)
    for key in CFG.realism.size_keys:
        assert bins[key]["bin_lo_bp"] == CFG.realism.bin_lo_bp()
        assert len(bins[key]["p"]) == CFG.realism.n_bins


def test_subset_of_known_names_is_accepted(tmp_path: Path) -> None:
    data = {
        "schema_version": 1,
        "design": {"compositions": {"markov": 1.0}, "delta_ranges": {"1": [1, 1], "2": [2, 2]}},
    }
    cfg = load_bench_config(_write(tmp_path, data))
    assert cfg.design.compositions == {"markov": 1.0} and list(cfg.design.delta_ranges) == [
        "1",
        "2",
    ]
