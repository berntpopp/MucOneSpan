"""Task 15e fix round 1: the single-event split needs peak-level evidence for its share.

After a single-event split the reads of each allele are selected by the event, so the
event's read support is circular. The split is therefore made only when the one-sided
lower confidence bound (``phase_single_event_alpha``) of the stutter-deconvolved minor
share reaches ``het_af_min``; otherwise the site stays unresolved (INCONCLUSIVE with its
located repeat), never PATHOGENIC.
"""

from __future__ import annotations

import dataclasses
import random
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.hybrid import single_event
from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.spans import ReadRecord, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_single_event as base

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
S = base.S
REPEAT = len(synth.PRE) + 1  # 1-based repeat of the first inner unit
# Realistic ONT amplicon depth for the stress cell (reads on the one length peak).
STRESS_DEPTH = 4 * S.phase_single_event_bound_reads
# Site-specific +1 excess at one C7 run, as a share of all reads, from well below to
# above het_af_min (both strands alike, so no strand test can reject it).
EXCESS_STEPS = (-2, -1, 0, 1, 2)
EXCESS = tuple(round(S.het_af_min + k * S.het_af_min / 4, 3) for k in EXCESS_STEPS)
STRESS_UNIT = 10  # 0-based inner unit whose C7 run carries the excess
# Stutter regimes: the ONT-like strand-asymmetric profile of the dupA tests, and a light
# symmetric one (its "-" strand on both strands), where the x4 run-background floor
# is low and only het_af_min and the share bound stand before a split.
STUTTER = {
    "ont_asymmetric": base.ONT_LIKE_STUTTER,
    "light": {"+": base.ONT_LIKE_STUTTER["-"], "-": base.ONT_LIKE_STUTTER["-"]},
}


def _stress_records(excess: float, stutter: str, seed: int) -> list[ReadRecord]:
    """Wild-type X reads with ONT-like stutter and a +1 C excess at one C7 run."""
    cons = synth.allele(["X"] * 30)
    start = base._unit_run(cons, STRESS_UNIT)
    edited = cons[:start] + "C" + cons[start:]
    rng = random.Random(seed)
    out = []
    for i in range(STRESS_DEPTH):
        strand = "+" if rng.random() < 1 / 2 else "-"
        source = edited if rng.random() < excess else cons
        template = base._stutter(source, rng, *STUTTER[stutter][strand])
        read = synth.reads(template, 1, err=base.ERR, seed=rng.randrange(1 << 30), strand_mix=False)
        seq = read[0].seq if strand == "+" else rc(read[0].seq)
        out.append(ReadRecord(f"x{seed}_{i}", seq, read[0].qual))
    return out


@pytest.mark.parametrize("stutter", sorted(STUTTER))
@pytest.mark.parametrize("excess", EXCESS)
def test_wild_type_site_specific_excess_is_never_pathogenic(
    tmp_path: Path, excess: float, stutter: str
) -> None:
    summary, decision = base._run(tmp_path, _stress_records(excess, stutter, S.seed))
    assert decision["state"] != "PATHOGENIC", (excess, summary["hybrid"]["split_bases"])
    if summary["hybrid"]["split_bases"] != ["none"]:
        assert decision["state"] == "INCONCLUSIVE"
        assert any("unresolved heterozygous site at repeat" in d for d in decision["details"])


def test_low_share_single_event_is_not_split() -> None:
    """A minor share just above het_af_min is not significantly above it: no split.

    A one-unit substitution (X -> Q, split only in "all" mode) has no stutter, so the
    observed share is the true share.
    """
    n = S.phase_max_site_reads
    minor = round(n * (S.het_af_min + S.het_af_min / 4))
    one_q = synth.allele(["X"] * 10 + ["Q"] + ["X"] * 19)
    wild = synth.allele(["X"] * 30)
    reads = synth.reads(one_q, minor, err=base.ERR, seed=3)
    reads += synth.reads(wild, n - minor, err=base.ERR, seed=4)
    members = categorize_reads(reads, base.ANCH, S).spanning
    s = dataclasses.replace(S, phase_single_event_split="all")
    res = split_by_linked_sites(wild, members, s, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site", res.basis
    assert single_event.split_single_event(wild, members, res, s) is None
    loose = dataclasses.replace(s, phase_single_event_alpha=1 / 2)  # bound = point estimate
    assert single_event.split_single_event(wild, members, res, loose) is not None


def test_single_event_alpha_is_validated() -> None:
    with pytest.raises(ValueError, match="phase_single_event_alpha"):
        dataclasses.replace(S, phase_single_event_alpha=1.0)


def test_unsplit_site_reason_reaches_the_decision(tmp_path: Path) -> None:
    settings = dataclasses.replace(
        DEFAULT_SETTINGS,
        hybrid=dataclasses.replace(DEFAULT_SETTINGS.hybrid, phase_single_event_split="off"),
    )
    _summary, decision = base._run(tmp_path, base._het_records(base.HET_SEEDS[0]), settings)
    assert decision["state"] == "INCONCLUSIVE"
    text = f"unresolved heterozygous site at repeat {REPEAT}"
    assert any(text in d for d in decision["details"]), decision["details"]


def test_equally_likely_run_lengths_leave_a_read_unassigned() -> None:
    """Observed 7 between major 6 and minor 8 under a symmetric error profile: a tie."""
    cap = S.phase_run_error_cap
    symmetric = [1 / (2 * cap + 1)] * (2 * cap + 1)
    profile = {"+": symmetric}
    site = {"site": ("run", 0), "major": 6, "minor": 8}
    assert single_event._allele({("run", 0): 7}, "+", site, (profile, profile)) is None
    assert single_event._allele({("run", 0): 8}, "+", site, None) == 1
    assert single_event._allele({("run", 0): 5}, "+", site, None) is None


def test_unassigned_reads_of_a_single_event_split_are_reassigned(tmp_path: Path) -> None:
    """Reads the split leaves unassigned (ties) are placed by edit distance, not dropped."""
    pytest.importorskip("pyabpoa", reason="pyabpoa (extra 'hybrid') is not installed")
    records = base._het_records(base.HET_SEEDS[0])
    tied = {r.name for r in records[:: len(records) // S.phase_min_minor_reads]}
    real = single_event.split_single_event

    def with_ties(cons: str, members: list[Any], res: Any, settings: Any) -> Any:
        out = real(cons, members, res, settings)
        if out is None:
            return None
        keep = [[m for m in g if m.name not in tied] for g in out.groups]
        moved = [m for g in out.groups for m in g if m.name in tied]
        return dataclasses.replace(out, groups=keep, unassigned=[*out.unassigned, *moved])

    fq = base._fastq(tmp_path / "in.fastq", records)
    with patch("muc_one_span.hybrid.engine.split_single_event", with_ties):
        result = reconstruct_alleles(fq, tmp_path, synth.RD, DEFAULT_SETTINGS)
    assert result.block["split_bases"] == ["single_event"]
    assert result.block["phase_unassigned_spanning_reads"] < len(tied)
    assigned = sum(result.alleles[k]["spanning_reads"] for k in ("allele_1", "allele_2"))
    total = result.block["read_categories"]["spanning"]
    assert assigned == total - result.block["unassigned_spanning_reads"]


def test_share_bound_sample_is_its_own_setting() -> None:
    """Raising the site-table compute cap must not weaken the share gate.

    At four times the bound's sample size a share of 1.25 x het_af_min is significant
    over every read but not over the bound's fixed sample.
    """
    n = 4 * S.phase_single_event_bound_reads
    minor = round(n * (S.het_af_min + S.het_af_min / 4))
    one_q = synth.allele(["X"] * 10 + ["Q"] + ["X"] * 19)
    wild = synth.allele(["X"] * 30)
    reads = synth.reads(one_q, minor, err=base.ERR, seed=5)
    reads += synth.reads(wild, n - minor, err=base.ERR, seed=6)
    members = categorize_reads(reads, base.ANCH, S).spanning
    s = dataclasses.replace(S, phase_single_event_split="all")
    res = split_by_linked_sites(wild, members, s, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site", res.basis
    assert single_event.split_single_event(wild, members, res, s) is None
    big_cap = dataclasses.replace(s, phase_max_site_reads=len(members))
    assert single_event.split_single_event(wild, members, res, big_cap) is None
    every_read = dataclasses.replace(s, phase_single_event_bound_reads=len(members))
    assert single_event.split_single_event(wild, members, res, every_read) is not None


def test_share_bound_reads_is_validated() -> None:
    with pytest.raises(ValueError, match="phase_single_event_bound_reads"):
        dataclasses.replace(S, phase_single_event_bound_reads=0)
