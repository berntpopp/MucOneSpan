"""Task 15e: an equal-length heterozygote with one indel event must never be NEGATIVE.

Reproduces the shape of a simulated ONT amplicon case: both alleles have the same
number of repeat units, one carries dupA in its first repeat (the dictionary X + A,
one base longer than the wild-type C unit, which differs from it only by the C-run
length), and the reads carry strand-asymmetric homopolymer stutter (heavy on "+"
reads, light on "-" reads). Before the fix the exact Fisher strand-bias test read the
stutter asymmetry as strand bias, dropped the only heterozygous site and merged both
haplotypes into one wild-type consensus (NEGATIVE).
"""

from __future__ import annotations

import dataclasses
import json
import random
import re
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.hybrid import engine
from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.phase import PhaseResult, split_by_linked_sites
from muc_one_span.hybrid.poa import get_backend
from muc_one_span.hybrid.polish import _runs
from muc_one_span.hybrid.single_event import is_indel, split_single_event
from muc_one_span.hybrid.spans import Anchors, ReadRecord, SpanRead, categorize_reads
from muc_one_span.pipeline import execute_pipeline
from muc_one_span.report import compute_clinical_decision
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S)
ERR = 0.02
N_PER_ALLELE = S.phase_max_site_reads // 2
DUPA = next(
    seq
    for seq, (parent, name) in synth.RD.mutated_sequences.items()
    if parent == "X" and name == "dupA"
)
WT = ["C"] + ["X"] * 29
MUT = [DUPA] + ["X"] * 29
# One-base stutter (deletion, insertion) probability per homopolymer run, by read
# strand: heavy on "+" reads and light on "-" reads, as measured at the C runs of the
# simulated ONT case (C7 read as C7 by ~53 % of "+" but ~91 % of "-" reads).
ONT_LIKE_STUTTER = {"+": (0.35, 0.10), "-": (0.03, 0.03)}
# Seeds whose heterozygote was dropped as strand-biased before the fix (4 of the first
# 10 in scratch sweeps), plus two that already passed.
HET_SEEDS = (0, 1, 5, 9, 2, 3)
WT_SEEDS = (0, 1, 2)
STRAND_ERROR_FRACS = (0.5, 0.7, 1.0)


def _stutter(seq: str, rng: random.Random, p_del: float, p_ins: float) -> str:
    out = []
    for match in re.finditer(r"(.)\1*", seq):
        run = match.group(0)
        if len(run) >= S.phase_run_min_len:
            r = rng.random()
            run = run[:-1] if r < p_del else run + run[0] if r < p_del + p_ins else run
        out.append(run)
    return "".join(out)


def _records(inner: list[str] | str, n: int, seed: int) -> list[ReadRecord]:
    """Stuttered reads; strand is drawn first so the stutter can depend on it."""
    rng = random.Random(seed)
    allele = synth.allele(inner) if isinstance(inner, list) else inner
    out = []
    for i in range(n):
        strand = "+" if rng.random() < 0.5 else "-"
        template = _stutter(allele, rng, *ONT_LIKE_STUTTER[strand])
        read = synth.reads(template, 1, err=ERR, seed=rng.randrange(1 << 30), strand_mix=False)[0]
        seq = read.seq if strand == "+" else rc(read.seq)
        out.append(ReadRecord(f"s{seed}_{i}", seq, read.qual))
    return out


def _spans(inner: list[str] | str, n: int, seed: int) -> list[SpanRead]:
    return categorize_reads(_records(inner, n, seed), ANCH, S).spanning


def _het(seed: int) -> list[SpanRead]:
    return _spans(MUT, N_PER_ALLELE, 2 * seed + 1) + _spans(WT, N_PER_ALLELE, 2 * seed + 2)


def _unit_run(cons: str, unit_index: int) -> int:
    """Start of the last C run of the given 0-based inner unit (after the pre-units).

    In both C and X that run is the one that dupA/dupC lengthen (C6 in C, C7 in X).
    """
    unit_bp = synth.RD.repeat_length_bp
    lo = (len(synth.PRE) + unit_index) * unit_bp
    runs = _runs(cons, S.phase_run_min_len)
    return max(s for s, _e, b in runs if b == "C" and lo <= s < lo + unit_bp)


@pytest.mark.parametrize("seed", HET_SEEDS)
def test_stutter_asymmetry_does_not_hide_a_heterozygous_run(seed: int) -> None:
    cons = synth.allele(WT)
    res = split_by_linked_sites(cons, _het(seed), S, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site", res.basis
    assert res.candidate is not None and res.candidate["site"] == ("run", _unit_run(cons, 0))


@pytest.mark.parametrize("seed", WT_SEEDS)
def test_wild_type_with_stutter_asymmetry_has_no_candidate(seed: int) -> None:
    cons = synth.allele(WT)
    members = _spans(WT, 2 * N_PER_ALLELE, 1000 + seed)
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    assert res.basis == "none"
    assert split_single_event(cons, members, res, S) is None


def _stranded_run_error(frac: float, seed: int) -> tuple[str, list[SpanRead]]:
    """Homozygote whose "-" reads carry C8 at one X unit's C7 run at fraction frac."""
    cons = synth.allele(["X"] * 30)
    start = _unit_run(cons, 10)
    edited = cons[:start] + "C" + cons[start:]
    rng = random.Random(seed)
    records = []
    for i in range(2 * N_PER_ALLELE):
        minus = rng.random() < 0.5
        source = edited if minus and rng.random() < frac else cons
        read = synth.reads(source, 1, err=ERR, seed=rng.randrange(1 << 30), strand_mix=False)[0]
        seq = rc(read.seq) if minus else read.seq
        records.append(ReadRecord(f"e{seed}_{i}", seq, read.qual))
    return cons, categorize_reads(records, ANCH, S).spanning


@pytest.mark.parametrize("frac", STRAND_ERROR_FRACS)
def test_strand_specific_run_error_is_still_rejected(frac: float) -> None:
    cons, members = _stranded_run_error(frac, S.seed)
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    assert res.basis == "none", res.sites


@pytest.mark.parametrize("seed", HET_SEEDS[:3])
def test_single_indel_event_splits_into_two_haplotypes(seed: int) -> None:
    cons = synth.allele(WT)
    first = _spans(MUT, N_PER_ALLELE, 2 * seed + 1)
    members = first + _spans(WT, N_PER_ALLELE, 2 * seed + 2)
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    split = split_single_event(cons, members, res, S)
    assert split is not None and split.basis == "single_event"
    names = {m.name for m in first}
    origin = [sum(m.name in names for m in g) / len(g) for g in split.groups]
    # Stutter moves some carrier reads to the other run length, so groups are not pure,
    # but each group is dominated by one haplotype and they are different ones.
    assert sorted(o > 1 / 2 for o in origin) == [False, True]
    assert min(len(g) for g in split.groups) >= S.het_min_group * len(members)
    grouped = [m.name for g in split.groups for m in g] + [m.name for m in split.unassigned]
    assert sorted(grouped) == sorted(m.name for m in members)


def _one_q() -> tuple[str, list[SpanRead]]:
    a = synth.allele(["X"] * 10 + ["Q"] + ["X"] * 19)
    reads = synth.reads(a, 40, err=ERR, seed=3) + synth.reads(
        synth.allele(["X"] * 30), 40, err=ERR, seed=4
    )
    return synth.allele(["X"] * 30), categorize_reads(reads, ANCH, S).spanning


@pytest.mark.parametrize(("mode", "splits"), [("off", False), ("indel", False), ("all", True)])
def test_single_event_mode_is_read_from_settings(mode: str, splits: bool) -> None:
    cons, members = _one_q()
    s = dataclasses.replace(S, phase_single_event_split=mode)
    res = split_by_linked_sites(cons, members, s, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site"
    assert (split_single_event(cons, members, res, s) is not None) is splits


def test_indel_split_is_off_when_disabled() -> None:
    cons, members = synth.allele(WT), _het(HET_SEEDS[0])
    s = dataclasses.replace(S, phase_single_event_split="off")
    res = split_by_linked_sites(cons, members, s, random.Random(S.seed))
    assert split_single_event(cons, members, res, s) is None


def test_single_event_mode_is_validated() -> None:
    with pytest.raises(ValueError, match="phase_single_event_split"):
        HybridSettings(phase_single_event_split="sometimes")


# --- end to end: the clinical decision ------------------------------------------------


def _run(tmp_path: Path, records: list[ReadRecord], settings: Any = DEFAULT_SETTINGS) -> Any:
    fq = tmp_path / "in.fastq"
    fq.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in records))
    out = tmp_path / "out"
    with patch("muc_one_span.tools.check_tools"):
        execute_pipeline(
            str(fq),
            str(out),
            None,
            "",
            1,
            10,
            5.0,
            False,
            "ont",
            None,
            engine="hybrid",
            settings=settings,
        )
    summary = json.loads((out / "summary.json").read_text())
    return summary, compute_clinical_decision(summary)


def _het_records(seed: int) -> list[ReadRecord]:
    return _records(MUT, N_PER_ALLELE, 2 * seed + 1) + _records(WT, N_PER_ALLELE, 2 * seed + 2)


@pytest.mark.parametrize("seed", HET_SEEDS[:2])
def test_equal_length_dupa_heterozygote_is_pathogenic(tmp_path: Path, seed: int) -> None:
    summary, decision = _run(tmp_path, _het_records(seed))
    assert summary["hybrid"]["split_bases"] == ["single_event"]
    mutations = [m for c in summary["classifications"].values() for m in c["mutations"]]
    dupa = [m for m in mutations if m.get("mutation_name") == "dupA"]
    assert len(dupa) == 1 and dupa[0]["read_support"]["status"] == "supported", dupa
    assert decision["state"] == "PATHOGENIC"


@pytest.mark.parametrize("seed", WT_SEEDS[:2])
def test_wild_type_with_realistic_stutter_stays_negative(tmp_path: Path, seed: int) -> None:
    summary, decision = _run(tmp_path, _records(WT, 2 * N_PER_ALLELE, 1000 + seed))
    assert summary["hybrid"]["split_bases"] == ["none"]
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"


def test_unsplit_heterozygous_site_is_inconclusive_with_its_repeat(tmp_path: Path) -> None:
    settings = dataclasses.replace(
        DEFAULT_SETTINGS,
        hybrid=dataclasses.replace(DEFAULT_SETTINGS.hybrid, phase_single_event_split="off"),
    )
    summary, decision = _run(tmp_path, _het_records(HET_SEEDS[0]), settings)
    assert summary["hybrid"]["split_bases"] == ["unconfirmed_single_site"]
    assert decision["state"] == "INCONCLUSIVE"
    repeat = len(synth.PRE) + 1  # the first inner unit, 1-based
    assert (
        f"unresolved heterozygous site at repeat {repeat}" in summary["hybrid"]["selection_detail"]
    )


# --- engine guards ----------------------------------------------------------------------


def _fastq(path: Path, records: list[ReadRecord]) -> Path:
    path.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in records))
    return path


def test_second_length_peak_is_not_split_on_a_single_event(tmp_path: Path) -> None:
    """With two length peaks each peak is one allele; a single event does not split it."""
    short = synth.reads(synth.allele(["X"] * 20), 2 * N_PER_ALLELE, err=ERR, seed=11)
    records = short + _het_records(HET_SEEDS[0])
    result = reconstruct_alleles(
        _fastq(tmp_path / "in.fastq", records), tmp_path, synth.RD, DEFAULT_SETTINGS
    )
    assert "single_event" not in result.block["split_bases"], result.block["split_bases"]
    assert "unconfirmed_single_site" in result.block["split_bases"]
    assert result.block["selection_status"] == "unresolved_single_site"
    repeat = len(synth.PRE) + 1
    assert f"unresolved heterozygous site at repeat {repeat}" in result.block["selection_detail"]


def test_identical_drafts_keep_the_peak_unconfirmed() -> None:
    """A single-event split whose two drafts are equal is not used (no haplotype split)."""
    members = _spans(WT, 2 * N_PER_ALLELE, 1000)
    cons = synth.allele(WT)
    unsplit = PhaseResult([members], "unconfirmed_single_site", [], candidate=None)
    half = len(members) // 2
    fake = PhaseResult([members[:half], members[half:]], "single_event")
    backend = get_backend(S.poa_backend)
    with patch("muc_one_span.hybrid.engine.split_single_event", return_value=fake):
        split, sub = engine._single_event(cons, members, unsplit, S, random.Random(S.seed), backend)
    assert sub is None and split is unsplit


def test_identical_polished_alleles_block_a_negative_call(tmp_path: Path) -> None:
    """If both single-event alleles polish to one sequence, the site stays unresolved."""
    real = engine.polish
    wild_type = synth.allele(WT)

    def same(draft: str, *args: Any, **kwargs: Any) -> tuple[str, dict[str, Any]]:
        return wild_type, real(draft, *args, **kwargs)[1]

    fq = _fastq(tmp_path / "in.fastq", _het_records(HET_SEEDS[0]))
    with patch("muc_one_span.hybrid.engine.polish", same):
        result = reconstruct_alleles(fq, tmp_path, synth.RD, DEFAULT_SETTINGS)
    assert result.block["split_bases"] == ["single_event"]
    assert result.block["selection_status"] == "unresolved_single_site"
    repeat = len(synth.PRE) + 1
    assert f"unresolved heterozygous site at repeat {repeat}" in result.block["selection_detail"]
    assert result.alleles["allele_1"]["selection_status"] == "unresolved_single_site"


def test_single_event_split_needs_both_groups_above_het_min_group() -> None:
    cons, members = synth.allele(WT), _het(HET_SEEDS[0])
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    assert split_single_event(cons, members, res, S) is not None
    # No split can give both groups more than half of the members.
    strict = dataclasses.replace(S, het_min_group=0.6)
    assert split_single_event(cons, members, res, strict) is None


@pytest.mark.parametrize(
    ("site", "major", "minor", "indel"),
    [
        (("run", 5), 6, 7, True),
        (("ins", 5), "", "A", True),
        (("ins", 5), "A", "C", False),
        (("col", 5), "A", "-", True),
        (("col", 5), "A", "C", False),
    ],
)
def test_is_indel_by_site_kind(site: Any, major: Any, minor: Any, indel: bool) -> None:
    assert is_indel({"site": site, "major": major, "minor": minor}) is indel
