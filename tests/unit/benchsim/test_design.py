from collections import Counter

from muc_one_span.benchsim.design import (
    DEPTHS,
    NORMAL_FRACTION,
    PROFILES,
    Design,
    build_split,
    derive_seed,
    event_bounds,
)

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
        assert sum(d.event is None for d in rows) / len(rows) >= NORMAL_FRACTION
        assert all(d.depth in DEPTHS[profile] for d in rows)


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
    # Unit 1 holds part of the forward primer site and the last five units are the
    # conserved 6-9 tail; events there break amplicon extraction or truth.
    assert event_bounds(20) == (5, 15)
    for split, n in (("dev", 300), ("val", 300)):
        for d in build_split(split, n, "salt", MUTS):
            for hap, repeat in d.targets:
                lo, hi = event_bounds(d.lengths[hap - 1])
                assert lo <= repeat <= hi, (d.design_id, d.lengths, d.targets)
