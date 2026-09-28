"""Task 15i: a strand-biased heterozygous insertion must never leave a NEGATIVE call.

Reproduces the shape of a simulated ONT amplicon case (a regression shape from the
unsealed test split, used non-blind): both alleles have the same repeat count and the
carrier differs by dupA, one A appended to an X unit (its C7 run followed by ``A`` then
the next unit's ``G``), so the event is an insertion slot, not a homopolymer run of
``phase_run_min_len`` bases. ONT "+" reads stutter heavily at the C7 run next to it; a
stuttered carrier read aligns its extra ``A`` into the run, so the insertion is seen in
far fewer "+" than "-" carrier reads. The column strand-bias test then refuses the
site, no candidate is left, both alleles merge into one wild-type consensus and, before
this task, the sample was NEGATIVE (the 15e stutter-aware test covers run sites only,
and the 15g safety tier tests run sites only).

The safety tier now also keeps an unsplit equal-length peak from a negative call when
a column or insertion site clears the candidate allele-fraction floor but is refused
only for strand bias: INCONCLUSIVE with the located reason, never an event.
"""

from __future__ import annotations

import random
import re
from pathlib import Path

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.engine import BIASED_SITE_BASIS, reconstruct_alleles
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.phase_sites import features, strand_biased_sites
from muc_one_span.hybrid.spans import ReadRecord, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_single_event as base

# Heavy synthetic safety sweep: its own CI job and make test-unit (never skipped).
pytestmark = pytest.mark.safety_sweep

S = DEFAULT_SETTINGS.hybrid
N_UNITS = 30
EVENT_UNIT = N_UNITS - 1  # last inner unit, followed by unit 6 (starts with G)
REPEAT = len(synth.PRE) + EVENT_UNIT + 1  # its 1-based repeat number
WILD_TYPE = ["X"] * N_UNITS
DUPA = ["X"] * EVENT_UNIT + [base.DUPA]
# Reads per allele: more than the phase sample, as at the case's depth.
N_PER_ALLELE = S.phase_max_site_reads
# Seeds whose het insertion the strand test refused before the fix (seed 0 passed it
# and became a single candidate site); the end-to-end tests use all of them.
SEEDS = (1, 2, 3)
ALL_SEEDS = (0, *SEEDS)
# One-base (deletion, insertion) stutter of the C7 runs by read strand, as measured at
# the event's C7 run in the case: read short in ~45% of "+" but ~5% of "-" reads. Only
# the C7 runs stutter here, so every other site stays clean and the test isolates the
# insertion next to a stuttering run.
C7 = 7
CASE_STUTTER = {"+": (0.45, 0.05), "-": (0.05, 0.03)}


def _stutter(seq: str, rng: random.Random, p_del: float, p_ins: float) -> str:
    out = []
    for match in re.finditer(r"(.)\1*", seq):
        run = match.group(0)
        if run[0] == "C" and len(run) >= C7:
            r = rng.random()
            run = run[:-1] if r < p_del else run + run[0] if r < p_del + p_ins else run
        out.append(run)
    return "".join(out)


def _records(inner: list[str], n: int, seed: int) -> list[ReadRecord]:
    """Stuttered reads; strand is drawn first so the stutter can depend on it."""
    rng = random.Random(seed)
    allele = synth.allele(inner)
    out = []
    for i in range(n):
        strand = "+" if rng.random() < 1 / 2 else "-"
        template = _stutter(allele, rng, *CASE_STUTTER[strand])
        read = synth.reads(template, 1, err=base.ERR, seed=rng.randrange(1 << 30), strand_mix=False)
        seq = read[0].seq if strand == "+" else rc(read[0].seq)
        out.append(ReadRecord(f"b{seed}_{i}", seq, read[0].qual))
    return out


def _het(seed: int) -> list[ReadRecord]:
    return _records(DUPA, N_PER_ALLELE, 2 * seed + 1) + _records(
        WILD_TYPE, N_PER_ALLELE, 2 * seed + 2
    )


def _wild(seed: int) -> list[ReadRecord]:
    return _records(WILD_TYPE, 2 * N_PER_ALLELE, 1000 + seed)


def _insertion_slot(cons: str) -> int:
    """Insertion slot after the event unit's last base (before unit 6)."""
    return (len(synth.PRE) + N_UNITS) * synth.RD.repeat_length_bp


@pytest.mark.parametrize("seed", SEEDS)
def test_case_shape_leaves_no_candidate_but_a_biased_site(seed: int) -> None:
    """The reproduction: the het insertion is refused as strand-biased, nothing else."""
    cons = synth.allele(WILD_TYPE)
    members = categorize_reads(_het(seed), base.ANCH, S).spanning
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    assert res.basis == "none", res.sites
    assert res.strand_biased, "the refused heterozygous insertion must be kept for the gate"
    top = res.strand_biased[0]
    assert top["site"] == ("ins", _insertion_slot(cons)), top
    assert top["minor"] == "A"


@pytest.mark.parametrize("seed", ALL_SEEDS)
def test_equal_length_strand_biased_dupa_is_never_negative(tmp_path: Path, seed: int) -> None:
    summary, decision = base._run(tmp_path, _het(seed))
    assert decision["state"] != "NO_PATHOGENIC_VARIANT_DETECTED", summary["hybrid"]
    if seed in SEEDS:  # the strand-biased shape: INCONCLUSIVE at the located site
        assert decision["state"] == "INCONCLUSIVE"
        text = f"unresolved heterozygous site at repeat {REPEAT}"
        assert any(text in d for d in decision["details"]), decision["details"]


@pytest.mark.parametrize("seed", ALL_SEEDS)
def test_wild_type_with_the_same_stutter_stays_negative(tmp_path: Path, seed: int) -> None:
    summary, decision = base._run(tmp_path, _wild(seed))
    assert summary["hybrid"]["split_bases"] == ["none"], summary["hybrid"]
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED", decision["details"]


def test_biased_site_tier_is_located_and_creates_no_event(tmp_path: Path) -> None:
    result = reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", _het(SEEDS[0])), tmp_path, synth.RD, DEFAULT_SETTINGS
    )
    assert result.block["split_bases"] == [BIASED_SITE_BASIS]
    assert result.block["selection_status"] == "unresolved_strand_biased_site"
    assert f"unresolved heterozygous site at repeat {REPEAT}" in result.block["selection_detail"]
    allele = result.alleles["allele_1"]
    assert allele["phase_status"] == "unresolved_strand_biased_site"
    assert allele["independent_haplotype_evidence"] is False
    # One group only: the tier never splits the peak into a carrier allele.
    assert result.alleles["allele_2"]["candidate_duplicate_of"] == "allele_1"
    assert result.alleles["homozygous"] is False


def test_strand_biased_sites_keep_only_het_level_refused_columns() -> None:
    """Sites below the candidate floor (het_af_min) are never kept; wild type has none."""
    cons = synth.allele(WILD_TYPE)
    for records, expect in ((_wild(SEEDS[0]), False), (_het(SEEDS[0]), True)):
        members = categorize_reads(records, base.ANCH, S).spanning
        sample = random.Random(S.seed).sample(members, S.phase_max_site_reads)
        feats, meta = features(cons, [m.seq for m in sample], S)
        sites = strand_biased_sites(feats, [m.strand for m in sample], meta, S)
        assert bool(sites) is expect, sites
        assert all(site["af"] >= S.het_af_min for site in sites)


def test_two_length_peaks_never_use_the_biased_site_tier(tmp_path: Path) -> None:
    """With two length peaks each peak is one allele; a within-peak site is not a haplotype."""
    short = synth.reads(
        synth.allele(["X"] * (N_UNITS // 2)), 2 * N_PER_ALLELE, err=base.ERR, seed=7
    )
    result = reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", short + _het(SEEDS[0])),
        tmp_path,
        synth.RD,
        DEFAULT_SETTINGS,
    )
    assert BIASED_SITE_BASIS not in result.block["split_bases"], result.block["split_bases"]
