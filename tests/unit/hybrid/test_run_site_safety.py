"""Task 15g: a heterozygous homopolymer run below the split floor must block NEGATIVE.

Reproduces the shape of a simulated HiFi amplicon case: both alleles have the same
repeat count and differ only at a dupC (one X unit's C7 run read as C8). HiFi-like
stutter reads a C7 run one base long in about 14% of reads, so the run-site candidate
floor (``phase_run_bg_multiplier`` x that background, about 0.56) sits above the
heterozygous C8 share (about 0.42). The site never becomes a candidate, both alleles
merge into one wild-type consensus and, before this task, the sample was NEGATIVE.

The safety tier tests every run of an unsplit equal-length peak against its
leave-one-out peer background at a lower, separately configured level; a run that
passes it makes the sample INCONCLUSIVE with the located reason and never creates an
event. Wild-type reads under realistic HiFi and ONT stutter must stay NEGATIVE.
"""

from __future__ import annotations

import dataclasses
import random
import re
from collections.abc import Callable
from pathlib import Path

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.engine import RUN_SITE_BASIS, reconstruct_alleles
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.phase_sites import features, run_excess_sites
from muc_one_span.hybrid.spans import ReadRecord, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_single_event as base
from tests.unit.hybrid import test_stutter_guard as guard

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
S = DEFAULT_SETTINGS.hybrid
N_UNITS = 30
# 0-based inner X unit carrying the dupC: near the end of the array, as in the case.
EVENT_UNIT = N_UNITS - 5
REPEAT = len(synth.PRE) + EVENT_UNIT + 1  # its 1-based repeat number
# Reads per allele: more than the phase sample, as at the case's 1000x depth.
N_PER_ALLELE = S.phase_max_site_reads
SEEDS = (0, 1, 2)
WILD_TYPE = ["X"] * N_UNITS
DUPC = ["X"] * EVENT_UNIT + [synth.dupc()] + ["X"] * (N_UNITS - EVENT_UNIT - 1)
# HiFi-like stutter measured on the case: C7 runs one base long in ~14% of reads and
# short in ~8%; the true C8 run short in ~20% (more deletion stutter at a longer run).
# Shorter runs stutter less, by a factor GROWTH per base.
C7 = 7
HIFI_INS = 0.14
HIFI_DEL = {C7: 0.08, C7 + 1: 0.20}
GROWTH = 2.5
# Light symmetric ONT-like "-" strand of the saturating ONT "+" shape.
ONT_MINUS_DEL = 0.03


def _hifi(length: int) -> tuple[float, float]:
    if length >= C7:
        return HIFI_DEL.get(length, max(HIFI_DEL.values())), HIFI_INS
    scale = GROWTH ** (length - C7)
    return HIFI_DEL[C7] * scale, HIFI_INS * scale


def _ont_plus(length: int) -> tuple[float, float]:
    return guard._ont_plus(length), base.ONT_LIKE_STUTTER["+"][1]


def _ont_minus(_length: int) -> tuple[float, float]:
    return ONT_MINUS_DEL, ONT_MINUS_DEL


Profile = Callable[[int], tuple[float, float]]
SHAPES: dict[str, dict[str, Profile]] = {
    "hifi": {"+": _hifi, "-": _hifi},
    "ont_saturating_plus": {"+": _ont_plus, "-": _ont_minus},
    "ont_asymmetric": {
        "+": lambda _n: base.ONT_LIKE_STUTTER["+"],
        "-": lambda _n: base.ONT_LIKE_STUTTER["-"],
    },
}


def _stutter(seq: str, rng: random.Random, profile: Profile) -> str:
    """Each run of >= phase_run_min_len bases loses or gains one base by its length."""
    out = []
    for match in re.finditer(r"(.)\1*", seq):
        run = match.group(0)
        if len(run) >= S.phase_run_min_len:
            p_del, p_ins = profile(len(run))
            r = rng.random()
            run = run[:-1] if r < p_del else run + run[0] if r < p_del + p_ins else run
        out.append(run)
    return "".join(out)


def _records(inner: list[str], n: int, seed: int, shape: str) -> list[ReadRecord]:
    rng = random.Random(seed)
    allele = synth.allele(inner)
    out = []
    for i in range(n):
        strand = "+" if rng.random() < 1 / 2 else "-"
        template = _stutter(allele, rng, SHAPES[shape][strand])
        read = synth.reads(template, 1, err=base.ERR, seed=rng.randrange(1 << 30), strand_mix=False)
        seq = read[0].seq if strand == "+" else rc(read[0].seq)
        out.append(ReadRecord(f"g{seed}_{i}", seq, read[0].qual))
    return out


def _het(seed: int, shape: str = "hifi") -> list[ReadRecord]:
    return _records(DUPC, N_PER_ALLELE, 2 * seed + 1, shape) + _records(
        WILD_TYPE, N_PER_ALLELE, 2 * seed + 2, shape
    )


def _wild(seed: int, shape: str) -> list[ReadRecord]:
    return _records(WILD_TYPE, 2 * N_PER_ALLELE, 1000 + seed, shape)


def test_case_shape_is_below_the_split_floor() -> None:
    """The reproduction: the het run is invisible to candidate-site detection."""
    members = categorize_reads(_het(SEEDS[0]), base.ANCH, S).spanning
    res = split_by_linked_sites(synth.allele(WILD_TYPE), members, S, random.Random(S.seed))
    assert res.basis == "none", res.sites


@pytest.mark.parametrize("seed", SEEDS)
def test_equal_length_dupc_below_the_split_floor_is_never_negative(
    tmp_path: Path, seed: int
) -> None:
    summary, decision = base._run(tmp_path, _het(seed))
    assert decision["state"] != "NO_PATHOGENIC_VARIANT_DETECTED", summary["hybrid"]
    if decision["state"] == "INCONCLUSIVE":
        text = f"unresolved heterozygous site at repeat {REPEAT}"
        assert any(text in d for d in decision["details"]), decision["details"]


@pytest.mark.parametrize("seed", SEEDS[:2])
@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_wild_type_under_realistic_stutter_stays_negative(
    tmp_path: Path, shape: str, seed: int
) -> None:
    summary, decision = base._run(tmp_path, _wild(seed, shape))
    assert summary["hybrid"]["split_bases"] == ["none"], summary["hybrid"]
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED", decision["details"]


def test_safety_tier_is_its_own_setting(tmp_path: Path) -> None:
    """At the split multiplier the tier adds nothing here: NEGATIVE, as before 15g."""
    off = dataclasses.replace(
        DEFAULT_SETTINGS,
        hybrid=dataclasses.replace(S, phase_run_safety_multiplier=S.phase_run_bg_multiplier),
    )
    _summary, decision = base._run(tmp_path, _het(SEEDS[0]), off)
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"


@pytest.mark.parametrize("seed", SEEDS[:1])
def test_run_site_tier_is_located_and_creates_no_event(tmp_path: Path, seed: int) -> None:
    result = reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", _het(seed)), tmp_path, synth.RD, DEFAULT_SETTINGS
    )
    assert result.block["split_bases"] == [RUN_SITE_BASIS]
    assert result.block["selection_status"] == "unresolved_run_site"
    assert f"unresolved heterozygous site at repeat {REPEAT}" in result.block["selection_detail"]
    allele = result.alleles["allele_1"]
    assert allele["phase_status"] == "unresolved_run_site"
    assert allele["independent_haplotype_evidence"] is False
    # One group only: the tier never splits the peak into a carrier allele.
    assert result.alleles["allele_2"]["candidate_duplicate_of"] == "allele_1"


def test_run_excess_sites_use_the_safety_multiplier() -> None:
    cons = synth.allele(WILD_TYPE)
    members = categorize_reads(_het(SEEDS[0]), base.ANCH, S).spanning
    sample = random.Random(S.seed).sample(members, S.phase_max_site_reads)
    feats, meta = features(cons, [m.seq for m in sample], S)
    sites = run_excess_sites(feats, meta, S)
    assert sites and sites[0]["site"] == ("run", base._unit_run(cons, EVENT_UNIT))
    assert (sites[0]["major"], sites[0]["minor"]) == (C7, C7 + 1)
    at_split_floor = dataclasses.replace(S, phase_run_safety_multiplier=S.phase_run_bg_multiplier)
    assert run_excess_sites(feats, meta, at_split_floor) == []


def test_two_length_peaks_never_use_the_run_site_tier(tmp_path: Path) -> None:
    """With two length peaks each peak is one allele; a sub-floor run is not a haplotype."""
    short = synth.reads(
        synth.allele(["X"] * (N_UNITS // 2)), 2 * N_PER_ALLELE, err=base.ERR, seed=7
    )
    result = reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", short + _het(SEEDS[0])),
        tmp_path,
        synth.RD,
        DEFAULT_SETTINGS,
    )
    assert RUN_SITE_BASIS not in result.block["split_bases"], result.block["split_bases"]
