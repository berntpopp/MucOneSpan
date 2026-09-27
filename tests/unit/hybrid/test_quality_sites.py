"""Task 15j: a low-accuracy read subset must not leave a phantom site in a length peak.

Reproduces the shape behind most HiFi normals reported INCONCLUSIVE on the development
panels: the length model finds two peaks (two alleles), and inside one peak a subset of
reads with lower base quality shares a systematic error at one site. The minor allele
reaches just above ``het_af_min`` (about 0.22), forms a single candidate site, and the
peak became ``unconfirmed_single_site`` although a within-peak site cannot be a further
allele when both alleles already have their own peak.

The rule (``phase_quality``) drops such a site from site detection only when its minor
carriers have significantly lower base quality than its major carriers
(``phase_quality_alpha``) and the minor allele fraction among the best
``phase_quality_keep_frac`` of reads is below ``het_af_min``. It applies only to peaks of
a two-peak model, never touches allele consensus, peak support or event evidence, and
fails closed when too few reads remain to re-test the site.
"""

from __future__ import annotations

import dataclasses
import math
import random
from pathlib import Path

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.phase_quality import (
    NO_QUALITY,
    binomial_lower_tail,
    quality_floor_reads,
    quality_sites,
    rank_sum_p_lower,
)
from muc_one_span.hybrid.phase_sites import candidates, features
from muc_one_span.hybrid.spans import PHRED_OFFSET, ReadRecord, SpanRead, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings, RuntimeSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_single_event as base

S = DEFAULT_SETTINGS.hybrid
SHORT = ["X"] * 20  # allele 1 (its own length peak)
LONG = ["X"] * 30  # allele 2; the artefact or minor haplotype sits inside its peak
ARTEFACT_UNIT = 15  # inner unit of LONG carrying the shared error / minor base
# Reads per allele: the phase sample size, as at the development panels' depth.
N_PER_ALLELE = S.phase_max_site_reads
# Share of the long allele's reads in the low-accuracy subset: its shared error lands
# at an allele fraction of about 0.22-0.25, just above het_af_min, as on the panels.
SUBSET_FRAC = 0.24
# Per-read mean Phred quality ranges (overlapping, as on the panels, where minor and
# major carriers differed by ~2 Q in their means). Low-accuracy reads also carry twice
# the random error rate of the others.
GOOD_Q = (30, 40)
POOR_Q = (22, 33)
POOR_ERR_FACTOR = 2
SEEDS = (0, 1)
# A permissive AF-bound level, to isolate the kept-read floor in its test.
LOOSE_AF_ALPHA = 0.5
# Units between the artefact unit and each of two linked artefact columns (non-touching
# events, so they count as two linked sites).
LINKED_UNIT_GAP = 5


def _column(template: str, unit: int = ARTEFACT_UNIT) -> tuple[int, str]:
    """A column of ``unit`` outside any run, and a base that makes no run there."""
    start = (len(synth.PRE) + unit) * synth.RD.repeat_length_bp
    for p in range(start + 1, start + synth.RD.repeat_length_bp - 1):
        near = {template[p - 1], template[p], template[p + 1]}
        if len(near) == 3:
            return p, next(b for b in "ACGT" if b not in near)
    raise AssertionError("no isolated column in the unit")


def _records(
    inner: list[str],
    n: int,
    seed: int,
    *,
    q: tuple[int, int],
    err: float,
    altered: bool,
    units: tuple[int, ...] = (ARTEFACT_UNIT,),
) -> list[ReadRecord]:
    """``n`` reads with per-read quality in ``q``; ``altered`` reads carry the site base
    in each of ``units``."""
    rng = random.Random(seed)
    template = synth.allele(inner)
    for unit in units if altered else ():
        pos, alt = _column(template, unit)
        template = template[:pos] + alt + template[pos + 1 :]
    out = []
    for i in range(n):
        read = synth.reads(template, 1, err=err, seed=rng.randrange(1 << 30), strand_mix=False)
        seq = read[0].seq if rng.random() < 1 / 2 else rc(read[0].seq)
        qual = chr(rng.randint(*q) + PHRED_OFFSET) * len(seq)
        out.append(ReadRecord(f"s{seed}_{int(altered)}_{i}", seq, qual))
    return out


def _sample(
    seed: int,
    *,
    n_long: int = N_PER_ALLELE,
    subset_q: tuple[int, int] = POOR_Q,
    subset_err_factor: int = POOR_ERR_FACTOR,
    long_inner: list[str] = LONG,
    poor_others: int = 0,
) -> list[ReadRecord]:
    """Two alleles; SUBSET_FRAC of the long allele's reads carry the site base.

    ``poor_others`` of the long allele's other reads are low-accuracy reads without it.
    """
    n_sub = round(SUBSET_FRAC * n_long)
    err = base.ERR
    n_good = n_long - n_sub - poor_others
    return (
        _records(SHORT, N_PER_ALLELE, 10 * seed + 1, q=GOOD_Q, err=err, altered=False)
        + _records(long_inner, n_good, 10 * seed + 2, q=GOOD_Q, err=err, altered=False)
        + _records(
            long_inner,
            poor_others,
            10 * seed + 4,
            q=POOR_Q,
            err=err * POOR_ERR_FACTOR,
            altered=False,
        )
        + _records(
            long_inner,
            n_sub,
            10 * seed + 3,
            q=subset_q,
            err=err * subset_err_factor,
            altered=True,
        )
    )


def _with(**changes: object) -> HybridSettings:
    return dataclasses.replace(S, **changes)  # type: ignore[arg-type]


def _long_peak(records: list[ReadRecord]) -> list[SpanRead]:
    spans = categorize_reads(records, base.ANCH, S).spanning
    cut = (len(synth.allele(SHORT)) + len(synth.allele(LONG))) // 2
    return [m for m in spans if m.length > cut]


def _settings_off() -> RuntimeSettings:
    return dataclasses.replace(DEFAULT_SETTINGS, hybrid=_with(phase_quality_alpha=0.0))


# --- the statistic ------------------------------------------------------------------


def test_rank_sum_p_is_small_only_when_the_first_sample_is_lower() -> None:
    low, high = [float(v) for v in range(10)], [float(v) for v in range(20, 40)]
    assert rank_sum_p_lower(low, high) < S.phase_quality_alpha
    assert rank_sum_p_lower(high, low) > 1 - S.phase_quality_alpha
    mixed = rank_sum_p_lower(low + high, high + low)
    assert S.phase_quality_alpha < mixed < 1 - S.phase_quality_alpha


def test_rank_sum_without_quality_variance_is_no_evidence() -> None:
    """Constant qualities (uninformative input) can never drop a site."""
    assert rank_sum_p_lower([30.0] * 20, [30.0] * 80) == 1.0
    assert rank_sum_p_lower([], [30.0, 31.0]) == 1.0


def test_binomial_lower_tail_is_exact() -> None:
    n, p = 10, S.het_af_min
    assert binomial_lower_tail(0, n, p) == pytest.approx((1 - p) ** n)
    assert binomial_lower_tail(n, n, p) == 1.0
    assert binomial_lower_tail(1, n, p) == pytest.approx((1 - p) ** n + n * p * (1 - p) ** (n - 1))


def test_a_read_without_qualities_keeps_every_site() -> None:
    """Mixed input: a quality-less read (mean Phred 0) would rank last; fail closed."""
    members = _long_peak(_sample(SEEDS[0]))
    cons = synth.allele(LONG)
    feats, meta = features(cons, [m.seq for m in members], S)
    sites = candidates(feats, [m.strand for m in members], meta, S)
    quals = [m.mean_q for m in members]
    assert quality_sites(sites, feats, quals, S)[1], "precondition: dropped with qualities"
    quals[0] = float(NO_QUALITY)
    kept, dropped = quality_sites(sites, feats, quals, S)
    assert kept == sites and not dropped


def test_quality_floor_is_derived_from_the_candidate_floors() -> None:
    """Enough kept reads to hold phase_min_minor_reads minor reads at het_af_min."""
    assert quality_floor_reads(S) == math.ceil(S.phase_min_minor_reads / S.het_af_min)


# --- (c) the artefact shape -----------------------------------------------------------


@pytest.mark.parametrize("seed", SEEDS)
def test_shape_is_a_candidate_without_the_rule(seed: int) -> None:
    """Precondition: the shared error is a candidate site of the long allele's peak."""
    members = _long_peak(_sample(seed))
    res = split_by_linked_sites(
        synth.allele(LONG), members, S, random.Random(S.seed), quality_filter=False
    )
    assert res.basis == "unconfirmed_single_site", res.sites
    assert S.het_af_min <= res.candidate["af"] < 2 * SUBSET_FRAC  # type: ignore[index]


@pytest.mark.parametrize("seed", SEEDS)
def test_low_accuracy_subset_site_is_dropped_in_a_two_peak_model(seed: int) -> None:
    members = _long_peak(_sample(seed))
    res = split_by_linked_sites(
        synth.allele(LONG), members, S, random.Random(S.seed), quality_filter=True
    )
    assert res.basis == "none", res.sites
    assert res.quality_associated, "the dropped site must be recorded"
    dropped = res.quality_associated[0]
    assert dropped["quality_p"] < S.phase_quality_alpha
    assert dropped["af_high_quality"] < S.het_af_min


@pytest.mark.parametrize("seed", SEEDS)
def test_low_accuracy_subset_artefact_leaves_a_resolved_negative(tmp_path: Path, seed: int) -> None:
    summary, decision = base._run(tmp_path, _sample(seed))
    block = summary["hybrid"]
    assert block["selection_status"] == "resolved", block["selection_detail"]
    assert block["split_bases"] == ["none", "none"]
    assert block["quality_associated_sites"], block
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED", decision["details"]


def test_rule_off_keeps_the_site_unresolved(tmp_path: Path) -> None:
    summary, decision = base._run(tmp_path, _sample(SEEDS[0]), _settings_off())
    assert summary["hybrid"]["selection_status"] == "unresolved_single_site"
    assert decision["state"] == "INCONCLUSIVE"


def test_a_linked_split_is_never_undone(tmp_path: Path) -> None:
    """The rule never changes group membership: a linked split stays, so every allele
    consensus is built from the same reads as without the rule (here the low-accuracy
    group stays a third, dropped group: unresolved_max_alleles)."""
    n_sub = round(SUBSET_FRAC * N_PER_ALLELE)
    linked_units = (ARTEFACT_UNIT - LINKED_UNIT_GAP, ARTEFACT_UNIT + LINKED_UNIT_GAP)
    records = (
        _records(SHORT, N_PER_ALLELE, 41, q=GOOD_Q, err=base.ERR, altered=False)
        + _records(LONG, N_PER_ALLELE - n_sub, 42, q=GOOD_Q, err=base.ERR, altered=False)
        + _records(
            LONG,
            n_sub,
            43,
            q=POOR_Q,
            err=POOR_ERR_FACTOR * base.ERR,
            altered=True,
            units=linked_units,
        )
    )
    members = _long_peak(records)
    res = split_by_linked_sites(
        synth.allele(LONG), members, S, random.Random(S.seed), quality_filter=True
    )
    assert res.basis == "linked_sites", res.sites
    assert not res.quality_associated
    # Task 15k: with phase_quality_group_exclusion (test_quality_groups) this poor group
    # is no longer counted as an allele; the split itself is unchanged either way.
    no_exclusion = dataclasses.replace(
        DEFAULT_SETTINGS, hybrid=_with(phase_quality_group_exclusion=False)
    )
    result = reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", records), tmp_path, synth.RD, no_exclusion
    )
    assert result.block["selection_status"] == "unresolved_max_alleles"


# --- (a) a true minor carried by good reads -------------------------------------------


@pytest.mark.parametrize("seed", SEEDS)
def test_minor_on_good_reads_is_still_a_site(tmp_path: Path, seed: int) -> None:
    records = _sample(seed, subset_q=GOOD_Q, subset_err_factor=1)
    summary, decision = base._run(tmp_path, records)
    assert summary["hybrid"]["selection_status"] == "unresolved_single_site"
    assert not summary["hybrid"]["quality_associated_sites"]
    assert decision["state"] == "INCONCLUSIVE"


@pytest.mark.parametrize("seed", SEEDS)
def test_minor_on_good_reads_survives_a_separate_poor_subset(seed: int) -> None:
    """Poor reads that do not carry the minor do not make its carriers look poor."""
    records = _sample(seed, subset_q=GOOD_Q, subset_err_factor=1, poor_others=N_PER_ALLELE // 4)
    members = _long_peak(records)
    res = split_by_linked_sites(
        synth.allele(LONG), members, S, random.Random(S.seed), quality_filter=True
    )
    assert res.basis == "unconfirmed_single_site", res.quality_associated
    assert not res.quality_associated


# --- (b) a true event on a real allele -------------------------------------------------


@pytest.mark.parametrize("seed", SEEDS)
def test_event_on_an_allele_is_still_called(tmp_path: Path, seed: int) -> None:
    """The rule never touches consensus or event evidence: a dupC allele stays called."""
    carrier = ["X"] * 5 + [synth.dupc()] + ["X"] * 24
    summary, decision = base._run(tmp_path, _sample(seed, long_inner=carrier))
    assert summary["hybrid"]["quality_associated_sites"], summary["hybrid"]
    mutations = [m for c in summary["classifications"].values() for m in c["mutations"]]
    assert any(m.get("mutation_name") == "dupC" for m in mutations), mutations
    assert decision["state"] == "PATHOGENIC", decision["details"]


def test_minor_event_on_good_reads_blocks_a_negative(tmp_path: Path) -> None:
    """A within-peak minor dupC carried by good reads is never released to NEGATIVE."""
    carrier = ["X"] * 5 + [synth.dupc()] + ["X"] * 24
    n_minor = round(SUBSET_FRAC * N_PER_ALLELE)
    n_poor = N_PER_ALLELE // 4
    records = (
        _records(SHORT, N_PER_ALLELE, 21, q=GOOD_Q, err=base.ERR, altered=False)
        + _records(LONG, N_PER_ALLELE - n_minor - n_poor, 22, q=GOOD_Q, err=base.ERR, altered=False)
        + _records(carrier, n_minor, 23, q=GOOD_Q, err=base.ERR, altered=False)
        + _records(LONG, n_poor, 24, q=POOR_Q, err=POOR_ERR_FACTOR * base.ERR, altered=False)
    )
    for sub in ("off", "on"):
        (tmp_path / sub).mkdir()
    off, off_decision = base._run(tmp_path / "off", records, _settings_off())
    assert off["hybrid"]["selection_status"] == "unresolved_single_site", "precondition"
    summary, decision = base._run(tmp_path / "on", records)
    assert not summary["hybrid"]["quality_associated_sites"]
    assert decision["state"] == off_decision["state"] == "INCONCLUSIVE", summary["hybrid"]


# --- (d) too few reads, uninformative qualities, scope --------------------------------


def test_too_few_kept_reads_fail_closed(tmp_path: Path) -> None:
    """Below the kept-read floor the site stays, even when the test would drop it."""
    # The largest peak whose high-quality subset (keep_frac of it) stays below the floor.
    n_long = math.floor((quality_floor_reads(S) - 1) / S.phase_quality_keep_frac)
    records = _sample(SEEDS[0], n_long=n_long)
    members = _long_peak(records)
    cons = synth.allele(LONG)
    feats, meta = features(cons, [m.seq for m in members], S)
    sites = candidates(feats, [m.strand for m in members], meta, S)
    assert sites, "precondition: the artefact is a candidate at this depth"
    quals = [m.mean_q for m in members]
    kept, dropped = quality_sites(sites, feats, quals, S)
    assert kept == sites and not dropped
    # The floor alone keeps it: with the AF bound made permissive the site is still
    # kept, and with a lower floor as well it is dropped.
    loose = _with(phase_quality_af_alpha=LOOSE_AF_ALPHA)
    assert not quality_sites(sites, feats, quals, loose)[1]
    lower = _with(phase_quality_af_alpha=LOOSE_AF_ALPHA, phase_min_minor_reads=1)
    assert quality_floor_reads(lower) < quality_floor_reads(S)
    assert quality_sites(sites, feats, quals, lower)[1]
    summary, decision = base._run(tmp_path, records)
    assert summary["hybrid"]["selection_status"] == "unresolved_single_site"
    assert decision["state"] == "INCONCLUSIVE"


def test_constant_qualities_keep_the_site() -> None:
    members = _long_peak(_sample(SEEDS[0], subset_q=GOOD_Q[:1] * 2))
    flat = [dataclasses.replace(m, mean_q=float(GOOD_Q[0])) for m in members]
    res = split_by_linked_sites(
        synth.allele(LONG), flat, S, random.Random(S.seed), quality_filter=True
    )
    assert res.basis == "unconfirmed_single_site"
    assert not res.quality_associated


def test_single_peak_is_out_of_scope(tmp_path: Path) -> None:
    """An equal-length genotype keeps every site: the rule is for two-peak models only."""
    n_sub = round(SUBSET_FRAC * 2 * N_PER_ALLELE)
    records = _records(
        LONG, 2 * N_PER_ALLELE - n_sub, 31, q=GOOD_Q, err=base.ERR, altered=False
    ) + _records(LONG, n_sub, 32, q=POOR_Q, err=POOR_ERR_FACTOR * base.ERR, altered=True)
    result = reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", records), tmp_path, synth.RD, DEFAULT_SETTINGS
    )
    assert result.block["split_bases"] == ["unconfirmed_single_site"]
    assert not result.block["quality_associated_sites"]


# --- settings ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("phase_quality_alpha", -0.1),
        ("phase_quality_alpha", 1.0),
        ("phase_quality_keep_frac", 0.0),
        ("phase_quality_af_alpha", 0.0),
        ("phase_quality_af_alpha", 1.0),
        ("phase_quality_keep_frac", 1.5),
    ],
)
def test_quality_settings_are_validated(field: str, value: float) -> None:
    with pytest.raises(ValueError, match=field):
        HybridSettings(**{field: value})  # type: ignore[arg-type]
