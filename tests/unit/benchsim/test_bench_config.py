import json
from dataclasses import replace
from pathlib import Path

import pytest

from muc_one_span.benchsim.bench_config import (
    DEFAULT_BENCH_CONFIG,
    BenchConfig,
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
    stress = json.loads(json.dumps(CFG.to_dict()))["sets"]["definitions"]["stress"]
    assert stress["offpeak_share_cap"] is None


def test_json_overlays_defaults_and_coerces_lists(tmp_path: Path) -> None:
    floor = CFG.amount.min_minor_share / 2
    path = _write(
        tmp_path,
        {
            "schema_version": 1,
            "amount": {"min_minor_share": floor},
            "design": {"length_min": CFG.design.length_min + 1},
        },
    )
    cfg = load_bench_config(path)
    assert cfg.amount.min_minor_share == floor
    assert cfg.design.length_min == CFG.design.length_min + 1
    assert cfg.realism == CFG.realism and cfg.sha256() != CFG.sha256()


@pytest.mark.parametrize(
    ("data", "match"),
    [
        ({"schema_version": 1, "nope": {}}, "unknown bench config sections"),
        ({"schema_version": 1, "report": {"nope": 1}}, "unknown report fields"),
        ({"report": {}}, "schema_version"),
        ({"schema_version": 1, "report": {"alpha": 0}}, "report.alpha"),
        ({"schema_version": 1, "structures": []}, "JSON object"),
        ({"schema_version": 1, "realism": {"offset_bin_min": 5}}, "offset_bin_min"),
        ({"schema_version": 1, "design": {"compositions": {"markov": 0.5}}}, "sum to 1"),
        ({"schema_version": 1, "amount": {"pcr_slope_per_unit": {"none": 0}}}, "pcr_slope"),
        ({"schema_version": 1, "design": {"delta_ranges": {">20": [21, 500]}}}, "length range"),
        ({"schema_version": 1, "design": {"split_sizes": {}}}, "at least one split"),
        ({"schema_version": 1, "design": {"depths": {}}}, "unknown design fields"),
        ({"schema_version": 1, "sets": {"definitions": []}}, "JSON object"),
        ({"schema_version": 1, "sets": {"definitions": {"standard": []}}}, "JSON object"),
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
        ({"schema_version": 1, "design": {"split_sizes": {"dev": "x"}}}, "split_sizes"),
        ({"schema_version": 1, "amount": {"pcr_slope_per_unit": {"none": "x"}}}, "pcr_slope"),
        ({"schema_version": 1, "atlas": {"decisions": ["MAYBE"]}}, r"atlas\.decisions"),
        ({"schema_version": 1, "atlas": {"decisions": []}}, r"atlas\.decisions"),
        ({"schema_version": 1, "atlas": {"strata": ["colour"]}}, r"atlas\.strata"),
        ({"schema_version": 1, "atlas": {"strata": [1]}}, r"atlas\.strata"),
        ({"schema_version": 1, "atlas": {"min_resolvable_depth": 0}}, "min_resolvable_depth"),
        ({"schema_version": 1, "atlas": {"min_resolvable_depth": 2.5}}, "min_resolvable_depth"),
        ({"schema_version": 1, "atlas": {"depth_basis": "median"}}, r"atlas\.depth_basis"),
        ({"schema_version": 1, "atlas": {"top_reasons": 0}}, r"atlas\.top_reasons"),
        (
            {"schema_version": 1, "atlas": {"expected_inconclusive_splits": ["nope"]}},
            "expected_inconclusive_splits",
        ),
        (
            {"schema_version": 1, "atlas": {"expected_inconclusive_splits": "stress"}},
            "JSON array",
        ),
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


def test_atlas_depth_gate_defaults_to_the_caller_per_allele_gate() -> None:
    from muc_one_span.settings import DEFAULT_SETTINGS

    gate = DEFAULT_SETTINGS.allele_selection.min_allele_primary_records
    assert CFG.atlas.min_resolvable_depth == gate
    assert "stress" in CFG.atlas.expected_inconclusive_splits
    assert set(CFG.atlas.expected_inconclusive_splits) <= set(CFG.design.split_sizes)


def test_atlas_overlay_accepts_known_names(tmp_path: Path) -> None:
    data = {
        "schema_version": 1,
        "atlas": {
            "decisions": ["INCONCLUSIVE", "NO_CALL"],
            "strata": ["depth"],
            "expected_inconclusive_splits": [],
            "depth_basis": "design",
            "min_resolvable_depth": CFG.atlas.min_resolvable_depth + 1,
        },
    }
    cfg = load_bench_config(_write(tmp_path, data))
    assert cfg.atlas.decisions == ("INCONCLUSIVE", "NO_CALL")
    assert cfg.atlas.expected_inconclusive_splits == ()
    assert cfg.atlas.depth_basis == "design" and cfg.sha256() != CFG.sha256()


def test_generation_hash_covers_only_generation_sections() -> None:
    from muc_one_span.benchsim.bench_config import GENERATION_SECTIONS

    assert set(GENERATION_SECTIONS) == {"design", "amount", "profiles", "structures"}
    same = replace(
        CFG,
        atlas=replace(CFG.atlas, top_reasons=CFG.atlas.top_reasons + 1),
        realism=replace(CFG.realism, jsd_max=CFG.realism.jsd_max / 2),
        run=replace(CFG.run, threads=CFG.run.threads + 1),
    )
    assert same.generation_sha256() == CFG.generation_sha256()
    assert same.sha256() != CFG.sha256()
    floor = CFG.amount.min_minor_share / 2
    changed = replace(CFG, amount=replace(CFG.amount, min_minor_share=floor))
    assert changed.generation_sha256() != CFG.generation_sha256()


def test_atlas_expected_sets_must_be_defined(tmp_path: Path) -> None:
    data = {"schema_version": 1, "atlas": {"expected_inconclusive_sets": ["nope"]}}
    with pytest.raises(ValueError, match="expected_inconclusive_sets"):
        load_bench_config(_write(tmp_path, data))
