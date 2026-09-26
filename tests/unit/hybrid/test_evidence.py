"""S8/S10 residual QC and read-level event support (homopolymer LLR and competition)."""

from __future__ import annotations

import random
from dataclasses import replace
from types import SimpleNamespace
from typing import Any

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.clinical_gates import READ_SUPPORT_STATUSES
from muc_one_span.hybrid.evidence import (
    competition_status,
    event_read_support,
    homopolymer_event_run,
    homopolymer_llr,
    hp_status,
    residual_sites,
)
from muc_one_span.hybrid.spans import Anchors, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
UNIT = synth.RD.repeat_length_bp
# Per-strand run-length profiles (index = observed length, 0..hp_max_run_len) around a
# 7-base reference run: the "+" strand stutters more than the "-" strand.
_TAIL = S.hp_max_run_len - 8
BACKGROUND = {
    "+": [0.02] * 6 + [0.20, 0.55, 0.20] + [0.03 / _TAIL] * _TAIL,
    "-": [0.01] * 6 + [0.05, 0.88, 0.05] + [0.01 / _TAIL] * _TAIL,
}
BOTH = {"+": 0, "-": 0}


def _oriented(allele: str, n: int, seed: int) -> list[tuple[str, str]]:
    cats = categorize_reads(
        synth.reads(allele, n, err=0.02, seed=seed), Anchors.from_dictionary(synth.RD, S), S
    )
    return [(sp.seq, sp.strand) for sp in cats.spanning]


def _variant(name: str, parent: str = "X") -> str:
    return next(
        seq for seq, (p, n) in synth.RD.mutated_sequences.items() if p == parent and n == name
    )


# --- residual QC ---------------------------------------------------------------------


def test_residual_sites_flag_mixture_but_not_clean_reads() -> None:
    cons = "ACGTACGTTTGACCATGCA" * 10
    alt = cons[:50] + "G" + cons[51:]
    af = S.qc_residual_af
    assert residual_sites(cons, [cons] * 20, af, min_run=S.qc_residual_min_run) == []
    sites = residual_sites(cons, [cons] * 12 + [alt] * 8, af, min_run=S.qc_residual_min_run)
    assert [(s["pos"], s["ref"], s["alt"], s["n"]) for s in sites] == [(50, "C", "G", 20)]


def test_residual_sites_skip_configured_runs() -> None:
    cons = "ACGTACGTTTGACCATGCA" * 10
    alt = cons[:8] + "C" + cons[9:]  # inside the TTT run at 7..10
    reads = [cons] * 12 + [alt] * 8
    assert residual_sites(cons, reads, S.qc_residual_af, min_run=S.qc_residual_min_run) == []
    flagged = residual_sites(cons, reads, S.qc_residual_af, min_run=S.qc_residual_min_run + 1)
    assert [s["pos"] for s in flagged] == [8]


# --- homopolymer model ---------------------------------------------------------------


def test_homopolymer_llr_separates_8c_from_7c_background() -> None:
    obs_mut = [("+", 8)] * 12 + [("+", 7)] * 8 + [("-", 8)] * 18 + [("-", 7)] * 2
    obs_wt = [("+", 7)] * 12 + [("+", 6)] * 8 + [("-", 7)] * 18 + [("-", 8)] * 2
    assert homopolymer_llr(obs_mut, BACKGROUND, 1) > S.hp_llr_min
    assert homopolymer_llr(obs_wt, BACKGROUND, 1) < 0


def test_homopolymer_llr_deletion_shift_mirrors_insertion() -> None:
    assert homopolymer_llr([("-", 6)] * 10, BACKGROUND, -1) > 0
    assert homopolymer_llr([("-", 8)] * 10, BACKGROUND, -1) < 0


def test_strand_without_background_contributes_nothing() -> None:
    assert homopolymer_llr([("-", 8)] * 5, {"+": BACKGROUND["+"]}, 1) == 0.0
    assert homopolymer_llr([], BACKGROUND, 1) == 0.0


# --- hp_status producer contract -------------------------------------------------------


def _ok_strands() -> tuple[dict[str, float], dict[str, int]]:
    n = S.hp_min_strand_reads
    return {"+": 1.0, "-": 1.0}, {"+": n, "-": n}


def test_hp_supported_when_every_threshold_is_met() -> None:
    llr_s, n_s = _ok_strands()
    assert hp_status(
        S.hp_min_reads, S.hp_llr_min, S.hp_min_alt_frac, llr_s, n_s, S, alternative_frac=0.0
    ) == ("supported")


def test_single_strand_input_is_supported_not_failed() -> None:
    obs = [("-", 8)] * S.hp_min_reads
    llr = homopolymer_llr(obs, BACKGROUND, 1)
    strand = {"+": homopolymer_llr([], BACKGROUND, 1), "-": llr}
    counts = {"+": 0, "-": S.hp_min_reads}
    assert (
        hp_status(S.hp_min_reads, llr, 1.0, strand, counts, S, alternative_frac=0.0) == "supported"
    )


def test_zero_read_strand_never_fails_even_with_negative_llr() -> None:
    llr_s, n_s = {"+": -5.0, "-": 3.0}, {"+": 0, "-": S.hp_min_reads}
    assert (
        hp_status(S.hp_min_reads, S.hp_llr_min, 1.0, llr_s, n_s, S, alternative_frac=0.0)
        == "supported"
    )


def test_thin_negative_strand_below_strand_minimum_does_not_fail() -> None:
    thin = S.hp_min_strand_reads - 1
    llr_s, n_s = {"+": -2.0, "-": 30.0}, {"+": thin, "-": S.hp_min_reads}
    assert (
        hp_status(S.hp_min_reads + thin, 30.0, 0.9, llr_s, n_s, S, alternative_frac=0.0)
        == "supported"
    )


def test_negative_well_covered_strand_is_discordant() -> None:
    plus = [("+", 7)] * S.hp_min_strand_reads
    minus = [("-", 8)] * S.hp_min_reads
    total = homopolymer_llr(plus + minus, BACKGROUND, 1)
    strand = {"+": homopolymer_llr(plus, BACKGROUND, 1), "-": homopolymer_llr(minus, BACKGROUND, 1)}
    counts = {"+": len(plus), "-": len(minus)}
    assert total > S.hp_llr_min and strand["+"] < 0
    assert (
        hp_status(len(plus + minus), total, 0.8, strand, counts, S, alternative_frac=0.0)
        == "discordant"
    )


def test_low_read_count_is_insufficient_depth() -> None:
    n = S.hp_min_reads - 1
    llr_s, n_s = _ok_strands()
    assert (
        hp_status(n, 10 * S.hp_llr_min, 1.0, llr_s, n_s, S, alternative_frac=0.0)
        == "insufficient_depth"
    )
    assert (
        hp_status(0, 0.0, 0.0, {"+": 0.0, "-": 0.0}, dict(BOTH), S, alternative_frac=0.0)
        == "insufficient_depth"
    )


@pytest.mark.parametrize("field", ["llr", "alt_frac"])
def test_hp_below_llr_or_alt_fraction_is_not_supported(field: str) -> None:
    llr_s, n_s = _ok_strands()
    llr = S.hp_llr_min - 0.01 if field == "llr" else S.hp_llr_min
    frac = S.hp_min_alt_frac - 0.01 if field == "alt_frac" else 1.0
    assert (
        hp_status(S.hp_min_reads, llr, frac, llr_s, n_s, S, alternative_frac=0.0) == "not_supported"
    )


def test_hp_thresholds_follow_settings() -> None:
    strict = replace(S, hp_min_reads=S.hp_min_reads + 5, hp_llr_min=S.hp_llr_min * 2)
    llr_s, n_s = _ok_strands()
    assert (
        hp_status(S.hp_min_reads, S.hp_llr_min, 1.0, llr_s, n_s, S, alternative_frac=0.0)
        == "supported"
    )
    assert hp_status(
        S.hp_min_reads, S.hp_llr_min, 1.0, llr_s, n_s, strict, alternative_frac=0.0
    ) == ("insufficient_depth")
    assert hp_status(
        strict.hp_min_reads, S.hp_llr_min, 1.0, llr_s, n_s, strict, alternative_frac=0.0
    ) == ("not_supported")


# --- competition producer contract ------------------------------------------------------


def test_competition_status_clauses() -> None:
    n = S.hp_min_reads
    alt = n  # all reads carry the template
    assert competition_status(n, alt, 0, S) == "supported"
    assert competition_status(n - 1, n - 1, 0, S) == "insufficient_depth"
    low_alt = int(S.hp_min_alt_frac * n) - 1
    assert competition_status(n, low_alt, 0, S) == "not_supported"
    half = n // 2
    assert competition_status(2 * half, half, half, S) == "not_supported"  # alt == ref


def _random_settings(rng: random.Random) -> HybridSettings:
    # Every other draw uses the most permissive valid configuration (boundaries).
    if rng.random() < 0.5:
        return replace(
            S,
            hp_min_reads=1,
            hp_llr_min=1e-9,
            hp_min_alt_frac=0.0,
            hp_min_strand_reads=0,
            event_max_alternative_frac=1.0,
        )
    return replace(
        S,
        hp_min_reads=rng.randint(1, 40),
        hp_llr_min=rng.uniform(1e-9, 30),
        hp_min_alt_frac=rng.uniform(0, 1),
        hp_min_strand_reads=rng.randint(0, 10),
        event_max_alternative_frac=rng.uniform(0, 1),
    )


def test_zero_reads_are_never_supported_even_at_permissive_settings() -> None:
    loose = replace(S, hp_min_reads=1, hp_llr_min=1e-9, hp_min_alt_frac=0.0)
    zero = {"+": 0, "-": 0}
    assert (
        hp_status(0, 50.0, 1.0, {"+": 1.0, "-": 1.0}, zero, loose, alternative_frac=0.0)
        == "insufficient_depth"
    )
    assert competition_status(0, 0, 0, loose) == "insufficient_depth"
    # The guard holds even for a settings object that bypassed validation.
    unchecked = SimpleNamespace(
        hp_min_reads=0, hp_llr_min=0.0, hp_min_alt_frac=0.0, hp_min_strand_reads=0
    )
    got = hp_status(0, 0.0, 0.0, {"+": 0.0, "-": 0.0}, zero, unchecked, alternative_frac=0.0)  # type: ignore[arg-type]
    assert got == "insufficient_depth"
    assert competition_status(0, 0, 0, unchecked) == "insufficient_depth"  # type: ignore[arg-type]


def test_property_status_set_and_thresholds() -> None:
    rng = random.Random(9)
    for _ in range(5000):
        s = _random_settings(rng)
        n = rng.randint(0, 60)
        llr = rng.uniform(-40, 60)
        frac = rng.uniform(0, 1)
        strand_n = {"+": rng.randint(0, n), "-": rng.randint(0, n)}
        strand_llr = {"+": rng.uniform(-20, 20), "-": rng.uniform(-20, 20)}
        mix = rng.uniform(0, 1)
        status = hp_status(n, llr, frac, strand_llr, strand_n, s, alternative_frac=mix)
        assert status in READ_SUPPORT_STATUSES
        assert (status == "insufficient_depth") == (n < max(s.hp_min_reads, 1))
        if status == "supported":
            assert (
                n >= 1 and n >= s.hp_min_reads and llr >= s.hp_llr_min and frac >= s.hp_min_alt_frac
            )
            assert mix <= s.event_max_alternative_frac
            assert not any(
                strand_n[st] >= s.hp_min_strand_reads and strand_n[st] > 0 and strand_llr[st] < 0
                for st in strand_n
            )
        alt = rng.randint(0, n)
        ref = rng.randint(0, n - alt)
        status = competition_status(n, alt, ref, s)
        assert status in READ_SUPPORT_STATUSES
        assert (status == "insufficient_depth") == (n < max(s.hp_min_reads, 1))
        if status == "supported":
            assert n >= 1 and n >= s.hp_min_reads and alt / n >= s.hp_min_alt_frac and alt > ref
            assert ref / n <= s.event_max_alternative_frac


# --- template typing -------------------------------------------------------------------


def test_event_type_comes_from_the_dictionary_template() -> None:
    cons = synth.allele(["X"] * 3 + [synth.dupc()] + ["X"] * 3)
    start = (len(synth.PRE) + 3) * UNIT
    end = start + UNIT + 1
    run = homopolymer_event_run({"mutation_name": "dupC"}, synth.RD, cons, start, end, S)
    x = synth.RD.repeats["X"]
    run_start = x.index("C" * 7)
    assert run == (start + run_start, start + run_start + 8, "C", 1)
    for name in ("dupA", "insG_pos58", "del18_31", "insCCCC", "delinsAT", "no_such"):
        mutation = {"mutation_name": name, "closest_type": "X"}
        assert homopolymer_event_run(mutation, synth.RD, cons, start, end, S) is None


def test_homopolymer_min_run_is_configured() -> None:
    cons = synth.allele(["X"] * 3 + [synth.dupc()] + ["X"] * 3)
    start = (len(synth.PRE) + 3) * UNIT
    event = {"mutation_name": "dupC"}
    too_long = replace(S, hp_event_min_run=9, hp_max_run_len=20)
    assert homopolymer_event_run(event, synth.RD, cons, start, start + UNIT + 1, too_long) is None
    capped = replace(S, hp_max_run_len=8)
    assert homopolymer_event_run(event, synth.RD, cons, start, start + UNIT + 1, capped) is None


def test_single_base_deletion_template_is_typed_from_the_parent() -> None:
    x = synth.RD.repeats["X"]
    run_start = x.index("C" * 7)
    template = {"changes": [{"type": "delete", "start": run_start + 3, "end": run_start + 3}]}
    rd = replace(synth.RD, mutations={**synth.RD.mutations, "delC_test": template})
    unit = x[: run_start + 2] + x[run_start + 3 :]
    cons = synth.allele(["X", unit, "X"])
    start = (len(synth.PRE) + 1) * UNIT
    mutation = {"mutation_name": "delC_test", "closest_type": "X"}
    got = homopolymer_event_run(mutation, rd, cons, start, start + UNIT - 1, S)
    assert got == (start + run_start, start + run_start + 6, "C", -1)
    assert homopolymer_event_run(dict(mutation, closest_type=None), rd, cons, 0, 1, S) is None


def test_run_crossing_the_unit_boundary_is_measured_whole() -> None:
    x = synth.RD.repeats["X"]
    parent = "TT" + x[2:]
    template = {"changes": [{"type": "insert", "start": 1, "sequence": "T"}]}
    rd = replace(synth.RD, mutations={**synth.RD.mutations, "insT_test": template})
    before = x[:-2] + "TT"  # the previous unit ends in the same base
    cons = synth.allele(["X", before, "T" + parent, "X"])
    start = (len(synth.PRE) + 2) * UNIT
    mutation = {"mutation_name": "insT_test", "closest_type": "X"}
    got = homopolymer_event_run(mutation, rd, cons, start, start + UNIT + 1, S)
    assert got == (start - 2, start + 3, "T", 1)


# --- end-to-end support -----------------------------------------------------------------


def _support(allele: str, reads: list[tuple[str, str]], name: str) -> dict[str, Any]:
    classification = classify_sequence(allele, synth.RD)
    idx = next(
        i for i, m in enumerate(classification["mutations_detected"]) if m["mutation_name"] == name
    )
    return event_read_support(allele, classification, reads, synth.RD, S)[idx]


def test_dupc_carrier_reads_support_the_event() -> None:
    allele = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    got = _support(allele, _oriented(allele, 40, 21), "dupC")
    assert got["kind"] == "homopolymer" and got["status"] == "supported", got
    assert got["n"] == got["alt"] + got["ref"] + got["other"]
    assert set(got["strand_llr"]) == {"+", "-"}


def test_wild_type_reads_do_not_support_a_consensus_dupc() -> None:
    allele = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    wild = _oriented(synth.allele(["X"] * 21), 40, 22)
    assert _support(allele, wild, "dupC")["status"] == "not_supported"


def test_non_homopolymer_event_uses_competition() -> None:
    allele = synth.allele(["X"] * 10 + [_variant("dupA")] + ["X"] * 10)
    got = _support(allele, _oriented(allele, 40, 23), "dupA")
    assert got["kind"] == "competition" and got["status"] == "supported", got


def test_wild_type_reads_do_not_support_a_competition_event() -> None:
    allele = synth.allele(["X"] * 10 + [_variant("dupA")] + ["X"] * 10)
    wild = _oriented(synth.allele(["X"] * 21), 40, 24)
    got = _support(allele, wild, "dupA")
    assert got["kind"] == "competition" and got["status"] == "not_supported", got


@pytest.mark.parametrize("name", ["dupC", "dupA"])
def test_few_reads_are_insufficient_depth_not_negative(name: str) -> None:
    unit = synth.dupc() if name == "dupC" else _variant(name)
    allele = synth.allele(["X"] * 10 + [unit] + ["X"] * 10)
    reads = _oriented(allele, 40, 25)[: S.hp_min_reads - 1]
    assert _support(allele, reads, name)["status"] == "insufficient_depth"


def test_unlocalized_event_is_reported_not_supported() -> None:
    allele = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    classification = classify_sequence(allele, synth.RD)
    for mut in classification["mutations_detected"]:
        mut["repeat_index"] = -1
    got = event_read_support(allele, classification, [], synth.RD, S)
    assert got and all(v["status"] == "not_localized" and v["n"] == 0 for v in got.values())
