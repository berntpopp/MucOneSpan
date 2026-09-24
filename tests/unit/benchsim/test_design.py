from collections import Counter
from dataclasses import replace

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, BenchConfig
from muc_one_span.benchsim.design import (
    CONSERVED_HEAD,
    CONSERVED_TAIL,
    PROFILES,
    Design,
    build_split,
    derive_seed,
    event_bounds,
)
from muc_one_span.settings import DEFAULT_LAYOUT

CFG = DEFAULT_BENCH_CONFIG.design
SETS = DEFAULT_BENCH_CONFIG.sets

MUTS = ["dupC", "dupA", "insG", "delGCCCA", "insCCC_benign"]


def test_seed_derivation_is_stable_and_stream_specific() -> None:
    a = derive_seed("s", "dev-0001", "bio")
    assert a == derive_seed("s", "dev-0001", "bio")
    assert a != derive_seed("s", "dev-0001", "reads") != derive_seed("t", "dev-0001", "bio")
    assert 0 <= a < 2**63


def test_split_sizes_normal_floor_and_profiles() -> None:
    designs = build_split("dev", 60, "salt", MUTS)
    assert len(designs) == 180 and {d.profile for d in designs} == set(PROFILES)
    for profile in PROFILES:
        rows = [d for d in designs if d.profile == profile]
        assert sum(d.event is None for d in rows) / len(rows) >= CFG.normal_fraction
        levels = SETS.definitions[SETS.default].profiles[profile]
        assert all(d.depth in levels.depths for d in rows)
        assert all(d.bench_set == SETS.default for d in rows)


def test_targets_match_event_allele_and_position() -> None:
    for d in build_split("dev", 90, "salt", MUTS):
        if d.event is None:
            assert d.targets == () and d.event_allele is None
            continue
        (hap, repeat), *_ = d.targets
        length = d.lengths[hap - 1]
        assert 1 <= repeat <= length
        if d.event_allele == "shorter" and d.lengths[0] != d.lengths[1]:
            assert length == min(d.lengths)
        if d.event_position == "first10":
            assert repeat <= max(1, length // 10) + 4  # MucOneUp allowed-parent search window


def test_designs_are_deterministic_and_round_trip() -> None:
    a, b = build_split("val", 30, "salt", MUTS), build_split("val", 30, "salt", MUTS)
    assert a == b
    assert Design.from_dict(a[0].to_dict()) == a[0]
    assert len({d.design_id for d in a}) == len(a)


def test_delta_classes_are_stratified() -> None:
    counts = Counter(d.delta_class for d in build_split("dev", 300, "salt", MUTS))
    assert min(counts.values()) >= 20


def test_event_targets_avoid_conserved_head_and_tail() -> None:
    # The conserved head (units 1-5, forward primer site in unit 1) and tail
    # (units 6-9) come from the bundled reference layout, shared with structures.
    assert (len(DEFAULT_LAYOUT.pre), len(DEFAULT_LAYOUT.after)) == (CONSERVED_HEAD, CONSERVED_TAIL)
    length = CONSERVED_HEAD + CONSERVED_TAIL + 1
    assert event_bounds(length) == (CONSERVED_HEAD + 1, CONSERVED_HEAD + 1)
    with pytest.raises(ValueError, match="too short"):
        event_bounds(length - 1)
    for split in ("dev", "val"):
        for d in build_split(split, CFG.split_sizes[split], "salt", MUTS):
            for hap, repeat in d.targets:
                lo, hi = event_bounds(d.lengths[hap - 1])
                assert lo <= repeat <= hi, (d.design_id, d.lengths, d.targets)


def test_clamped_targets_are_recorded() -> None:
    designs = build_split("dev", CFG.split_sizes["dev"], "salt", MUTS)
    clamped = [d for d in designs if d.target_clamped]
    assert clamped and all(d.event for d in clamped)
    assert not any(d.target_clamped for d in designs if d.event is None)
    for d in clamped:
        (hap, repeat), *_ = d.targets
        assert repeat in event_bounds(d.lengths[hap - 1])
    assert Design.from_dict(clamped[0].to_dict()) == clamped[0]
    legacy = {k: v for k, v in clamped[0].to_dict().items() if k != "target_clamped"}
    assert Design.from_dict(legacy).target_clamped is False


def test_factor_levels_follow_the_set() -> None:
    for name, bench_set in SETS.definitions.items():
        designs = build_split("dev", 30, "salt", MUTS, bench_set=name)
        for profile in PROFILES:
            rows = [d for d in designs if d.profile == profile]
            levels = bench_set.profiles[profile]
            assert {d.depth for d in rows} == set(levels.depths)
            assert {d.pcr for d in rows} == set(levels.pcr_levels)
            assert {d.error for d in rows} == set(levels.error_levels)
            assert {d.smear for d in rows} == set(levels.smear_levels)
            assert {d.chimera for d in rows} == set(levels.chimera_levels)
            assert all(d.bench_set == name for d in rows)
            assert all(d.design_id.startswith(f"dev-{name}-{profile}-") for d in rows)


def test_sets_share_the_biology_of_a_split() -> None:
    bio = ("profile", "lengths", "delta_class", "composition", "event", "event_allele")
    bio += ("event_position", "targets", "bio_seed", "target_clamped")
    runs = {name: build_split("dev", 30, "salt", MUTS, bench_set=name) for name in SETS.definitions}
    standard = runs.pop(SETS.default)
    for designs in runs.values():
        for a, b in zip(standard, designs, strict=True):
            assert [getattr(a, k) for k in bio] == [getattr(b, k) for k in bio]
            assert a.read_seed != b.read_seed and a.design_id != b.design_id


def test_unknown_set_is_rejected() -> None:
    with pytest.raises(ValueError, match="set"):
        build_split("dev", 3, "salt", MUTS, bench_set="nope")


def test_legacy_design_without_set_round_trips() -> None:
    design = build_split("dev", 1, "salt", MUTS)[0]
    legacy = {k: v for k, v in design.to_dict().items() if k != "bench_set"}
    assert Design.from_dict(legacy).bench_set is None


def test_design_factors_come_from_the_config() -> None:
    standard = SETS.definitions["standard"]
    profiles = {
        p: replace(levels, chimera_levels=(0.02,), pcr_levels=("none",))
        for p, levels in standard.profiles.items()
    }
    only = {"standard": replace(standard, profiles=profiles)}
    sets = replace(SETS, definitions=only, legacy="standard")
    cfg = BenchConfig(design=replace(CFG, normal_fraction=0.5, length_min=40), sets=sets)
    designs = build_split("val", 20, "salt", MUTS, cfg)
    assert {d.chimera for d in designs} == {0.02} and {d.pcr for d in designs} == {"none"}
    assert min(min(d.lengths) for d in designs) >= 40
    for profile in PROFILES:
        rows = [d for d in designs if d.profile == profile]
        assert sum(d.event is None for d in rows) == 10


def test_unknown_split_or_profile_levels_are_rejected() -> None:
    with pytest.raises(ValueError, match="split"):
        build_split("nope", 3, "salt", MUTS)
    standard = SETS.definitions["standard"]
    profiles = {k: v for k, v in standard.profiles.items() if k != "hifi_amplicon"}
    definitions = {**SETS.definitions, "standard": replace(standard, profiles=profiles)}
    cfg = BenchConfig(sets=replace(SETS, definitions=definitions))
    with pytest.raises(ValueError, match="hifi_amplicon"):
        build_split("dev", 3, "salt", MUTS, cfg)
