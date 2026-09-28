"""Task 15k: a linked group of low-accuracy reads is not a further allele.

Reproduces the validation HiFi normals reported INCONCLUSIVE (``unresolved_max_alleles``):
two length peaks, and inside one peak a subset of lower-quality reads shares two
systematic substitutions a few bases apart. The two sites are linked through the same
reads, so the peak splits into two groups and the sample holds three allele groups.
Task 15j keeps that split (the poor reads must never rejoin the allele consensus).

With ``phase_quality_group_exclusion`` the poor group still never joins an allele; it
only stops counting as a further allele when it is the smaller group, its reads have
lower base quality, its share of the highest-quality reads is significantly below
``het_af_min``, it has no distinct length and its draft differs from the other group's
by substitutions only (no event). A group carrying dupC, or carried by good reads,
keeps the sample unresolved.
"""

from __future__ import annotations

import dataclasses
import random
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.phase_groups import explained_group, substitutions_only
from muc_one_span.hybrid.phase_quality import NO_QUALITY, quality_floor_reads
from muc_one_span.hybrid.spans import ReadRecord, SpanRead
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings, RuntimeSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_quality_single_event as single
from tests.unit.hybrid import test_quality_sites as quality
from tests.unit.hybrid import test_single_event as base

# Both the 15j test and this rule are opt-in since Task 15k: switched on explicitly.
S = dataclasses.replace(quality.S, phase_quality_group_exclusion=True)
OPT_IN = dataclasses.replace(DEFAULT_SETTINGS, hybrid=S)
UNIT_BP = synth.RD.repeat_length_bp
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
N_PER_ALLELE = quality.N_PER_ALLELE
SEEDS = (0, 1)
LINKED_UNITS = (
    quality.ARTEFACT_UNIT - quality.LINKED_UNIT_GAP,
    quality.ARTEFACT_UNIT + quality.LINKED_UNIT_GAP,
)
# Minority group (share of the long allele's reads) for the adversarial sweep; its
# reads are low quality and carry the two linked substitutions, and with ``dupc`` also
# dupC at DUPC_UNIT. Low-quality non-carrier reads stutter C7 -> C8 at that run.
MINOR_AFS = (0.15, 0.20, 0.25, 0.30)
SWEEP_SEEDS = (0, 1, 2)
DUPC_UNIT = 3
POOR_STUTTER_AT_EVENT = single.POOR_STUTTER_AT_EVENT


def _settings(**changes: object) -> RuntimeSettings:
    hybrid: HybridSettings = dataclasses.replace(S, **changes)  # type: ignore[arg-type]
    return dataclasses.replace(DEFAULT_SETTINGS, hybrid=hybrid)


def _records(seed: int, *, subset_q: tuple[int, int] = quality.POOR_Q) -> list[ReadRecord]:
    """Two alleles; SUBSET_FRAC of the long allele's reads carry two linked columns."""
    n_sub = round(quality.SUBSET_FRAC * N_PER_ALLELE)
    return (
        quality._records(
            quality.SHORT,
            N_PER_ALLELE,
            10 * seed + 1,
            q=quality.GOOD_Q,
            err=base.ERR,
            altered=False,
        )
        + quality._records(
            quality.LONG,
            N_PER_ALLELE - n_sub,
            10 * seed + 2,
            q=quality.GOOD_Q,
            err=base.ERR,
            altered=False,
        )
        + quality._records(
            quality.LONG,
            n_sub,
            10 * seed + 3,
            q=subset_q,
            err=quality.POOR_ERR_FACTOR * base.ERR,
            altered=True,
            units=LINKED_UNITS,
        )
    )


def _reconstruct(tmp_path: Path, records: list[ReadRecord], settings: RuntimeSettings) -> Any:
    tmp_path.mkdir(exist_ok=True)
    return reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", records), tmp_path, synth.RD, settings
    )


# --- the rule, end to end -----------------------------------------------------------


@pytest.mark.parametrize("seed", SEEDS)
def test_poor_linked_group_is_not_a_third_allele(tmp_path: Path, seed: int) -> None:
    records = _records(seed)
    off = _reconstruct(tmp_path / "off", records, _settings(phase_quality_group_exclusion=False))
    assert off.block["selection_status"] == "unresolved_max_alleles", "precondition"
    on = _reconstruct(tmp_path / "on", records, OPT_IN)
    assert on.block["selection_status"] == "resolved", on.block["selection_detail"]
    excluded = on.block["quality_excluded_groups"]
    assert len(excluded) == 1 and excluded[0]["share_high_quality"] < S.het_af_min
    assert not on.block["dropped_groups"]
    # The poor reads never join an allele: both consensus sequences are unchanged.
    for name in ("allele_1", "allele_2"):
        assert on.consensus_paths[name].read_text() == off.consensus_paths[name].read_text()
    assert on.block["unassigned_spanning_reads"] >= excluded[0]["spanning_reads"]


def test_excluded_group_gives_a_negative_decision(tmp_path: Path) -> None:
    summary, decision = base._run(tmp_path, _records(SEEDS[0]), OPT_IN)
    assert summary["hybrid"]["quality_excluded_groups"], summary["hybrid"]
    assert decision["state"] == NEGATIVE, decision["details"]


def test_linked_group_on_good_reads_stays_a_third_allele(tmp_path: Path) -> None:
    result = _reconstruct(tmp_path, _records(SEEDS[0], subset_q=quality.GOOD_Q), OPT_IN)
    assert result.block["selection_status"] == "unresolved_max_alleles"
    assert not result.block["quality_excluded_groups"]


def test_rule_needs_the_quality_test() -> None:
    """Without the 15j test (phase_quality_alpha 0) the rule would be a silent no-op:
    the configuration is refused, naming both keys."""
    with pytest.raises(ValueError, match=r"phase_quality_group_exclusion.*phase_quality_alpha"):
        _settings(phase_quality_alpha=0.0)


def test_setting_is_validated() -> None:
    with pytest.raises(ValueError, match="phase_quality_group_exclusion"):
        HybridSettings(phase_quality_group_exclusion="yes")  # type: ignore[arg-type]


# --- the conditions, on the group test itself -----------------------------------------


def _span(i: int, seq: str, q: float) -> SpanRead:
    return SpanRead(f"g{i}", seq, q, "+", 0, "motif")


def _groups(
    small_seq: str, large_seq: str, *, small_q: float, n_small: int, n_large: int
) -> list[list[SpanRead]]:
    """A small group at quality ``small_q`` and a large one spread over GOOD_Q."""
    lo, hi = quality.GOOD_Q
    large = [_span(i, large_seq, lo + (hi - lo) * i / n_large) for i in range(n_large)]
    small = [_span(n_large + i, small_seq, small_q + i / n_small) for i in range(n_small)]
    return [large, small]


LARGE_SEQ = synth.allele(quality.LONG)
SUB_SEQ = quality._records(quality.LONG, 1, 0, q=quality.GOOD_Q, err=0.0, altered=True)[0].seq
N_LARGE = 3 * quality_floor_reads(S)
N_SMALL = N_LARGE // 3


def test_a_split_without_two_groups_is_kept() -> None:
    small_seq = LARGE_SEQ[:100] + ("A" if LARGE_SEQ[100] != "A" else "C") + LARGE_SEQ[101:]
    groups = _groups(
        small_seq, LARGE_SEQ, small_q=quality.POOR_Q[0], n_small=N_SMALL, n_large=N_LARGE
    )
    assert explained_group(groups[:1], [LARGE_SEQ], S, UNIT_BP) is None


def test_substitutions_only_needs_equal_length_and_no_indel() -> None:
    assert substitutions_only("ACGT", "ACCT")
    assert not substitutions_only("ACGT", "ACGGT")
    assert not substitutions_only("ACGTA", "CGTAA")  # shift: two indels beat four mismatches


def test_group_explained_by_poor_reads() -> None:
    small_seq = LARGE_SEQ[:100] + ("A" if LARGE_SEQ[100] != "A" else "C") + LARGE_SEQ[101:]
    groups = _groups(
        small_seq, LARGE_SEQ, small_q=quality.POOR_Q[0], n_small=N_SMALL, n_large=N_LARGE
    )
    hit = explained_group(groups, [LARGE_SEQ, small_seq], S, UNIT_BP)
    assert hit is not None and hit[0] == 1, hit


@pytest.mark.parametrize(
    "case", ["indel", "distinct_length", "no_quality", "good_quality", "too_few", "off"]
)
def test_group_is_kept_unless_every_condition_holds(case: str) -> None:
    small_seq = LARGE_SEQ[:100] + ("A" if LARGE_SEQ[100] != "A" else "C") + LARGE_SEQ[101:]
    small_q, n_small, n_large, s = float(quality.POOR_Q[0]), N_SMALL, N_LARGE, S
    if case == "indel":
        small_seq = LARGE_SEQ[:100] + LARGE_SEQ[100] + LARGE_SEQ[100:]  # a dupC-like insertion
    elif case == "distinct_length":
        small_seq = small_seq + LARGE_SEQ[-round(S.peak_min_separation_units * UNIT_BP) :]
    elif case == "no_quality":
        small_q = NO_QUALITY
    elif case == "good_quality":
        small_q = float(quality.GOOD_Q[1])
    elif case == "too_few":
        n_large = n_small = quality_floor_reads(S) // 2
    else:
        s = dataclasses.replace(S, phase_quality_group_exclusion=False)
    groups = _groups(small_seq, LARGE_SEQ, small_q=small_q, n_small=n_small, n_large=n_large)
    assert explained_group(groups, [LARGE_SEQ, small_seq], s, UNIT_BP) is None


# --- adversarial: a minority group carrying dupC --------------------------------------


def _minority(profile: str, af: float, seed: int, *, dupc: bool) -> list[ReadRecord]:
    """Two alleles; ``af`` of the long allele's reads form a low-quality minority with
    the two linked substitutions (and dupC at DUPC_UNIT when ``dupc``)."""
    rng = random.Random(seed)
    q = single.PROFILES[profile]
    inner = list(quality.LONG)
    if dupc:
        inner[DUPC_UNIT] = synth.dupc()
    minority = synth.allele(inner)
    for unit in LINKED_UNITS:
        pos, alt = quality._column(synth.allele(quality.LONG), unit)
        pos += len(minority) - len(synth.allele(quality.LONG))
        minority = minority[:pos] + alt + minority[pos + 1 :]
    stuttered = list(quality.LONG)
    stuttered[DUPC_UNIT] = synth.dupc()
    long_stutter = synth.allele(stuttered)
    out = []
    for i in range(N_PER_ALLELE):
        poor = rng.random() < single.POOR_FRAC
        out.append(
            quality_read(synth.allele(quality.SHORT), rng, f"a{i}", q["poor" if poor else "good"])
        )
    n_minor = round(af * N_PER_ALLELE)
    for i in range(N_PER_ALLELE):
        poor = i < n_minor or rng.random() < single.POOR_FRAC
        if i < n_minor:
            template = minority
        elif poor and dupc and rng.random() < POOR_STUTTER_AT_EVENT:
            template = long_stutter
        else:
            template = synth.allele(quality.LONG)
        out.append(quality_read(template, rng, f"b{i}", q["poor" if poor else "good"]))
    return out


def quality_read(template: str, rng: random.Random, name: str, q: tuple[int, int]) -> ReadRecord:
    strand = "+" if rng.random() < 1 / 2 else "-"
    return single._read(template, strand, rng.randrange(1 << 30), name, q)


# Task 15l: these shapes (two-peak dupC minority, all carriers low quality) were
# NEGATIVE-on-pathogenic strict xfails for the within-peak floor; the run-minority tier
# (hybrid.run_minor) now blocks them, so every case must pass.
DUPC_SWEEP = [
    (p, af, seed) for p in sorted(single.PROFILES) for af in MINOR_AFS for seed in SWEEP_SEEDS
]


@pytest.mark.parametrize(("profile", "af", "seed"), DUPC_SWEEP)
def test_dupc_minority_group_is_never_excluded(
    tmp_path: Path, profile: str, af: float, seed: int
) -> None:
    """A group carrying dupC is never excluded, and the sample is never NEGATIVE."""
    summary, decision = base._run(tmp_path, _minority(profile, af, seed, dupc=True), OPT_IN)
    assert not summary["hybrid"]["quality_excluded_groups"], summary["hybrid"]
    assert decision["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


@pytest.mark.parametrize("profile", sorted(single.PROFILES))
@pytest.mark.parametrize("af", MINOR_AFS)
@pytest.mark.parametrize("seed", SWEEP_SEEDS)
def test_wild_type_minority_group_is_never_pathogenic(
    tmp_path: Path, profile: str, af: float, seed: int
) -> None:
    summary, decision = base._run(tmp_path, _minority(profile, af, seed, dupc=False), OPT_IN)
    assert decision["state"] != "PATHOGENIC", summary["hybrid"]["split_bases"]
