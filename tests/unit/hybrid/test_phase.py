"""S4 linked-site phase split for equal and close allele lengths."""

from __future__ import annotations

import dataclasses
import math
import random
from collections.abc import Callable

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.phase import PhaseResult, split_by_linked_sites
from muc_one_span.hybrid.polish import _runs
from muc_one_span.hybrid.spans import Anchors, ReadRecord, SpanRead, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid.synth import LAYOUT

# Every phase tunable comes from S (HybridSettings defaults); nothing is duplicated here
# as a bare literal (owner directive: fully config-driven, no magic numbers).
S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S, LAYOUT)
N_PER_ALLELE = 40
ERR = 0.02
PLAIN = ["X"] * 30
# A and B differ from X at 4 and 2 non-adjacent bases; Q differs from X at one base.
TWO_UNITS = ["X"] * 10 + ["A"] + ["X"] * 10 + ["B"] + ["X"] * 8
ONE_A = ["X"] * 10 + ["A"] + ["X"] * 19
ONE_B = ["X"] * 10 + ["B"] + ["X"] * 19
ONE_Q = ["X"] * 10 + ["Q"] + ["X"] * 19
# Seed pairs (allele, PLAIN) sweeping the orientation cases the review found; the ONE_B
# pairs 302/303, 304/305 and 314/315 (in this sweep) gave opposite-oriented sites before
# fix round 1.
PURITY_SEEDS = [(300 + 2 * k, 301 + 2 * k) for k in range(8)]
# Fractions of "-" strand reads carrying a strand-specific systematic error (1.0: all).
STRAND_ERROR_FRACS = (0.5, 0.7, 0.9, 1.0)
# Unbiased imbalanced heterozygotes (variant : plain read ratios) must still split; fix
# round 1's per-strand AF floor dropped them to "none" by binomial strand sampling alone.
ALLELE_RATIOS = ((3, 1), (1, 3))
IMBALANCE_SEEDS = [(9500 + 2 * k, 9501 + 2 * k) for k in range(3)]


def _spans(inner: list[str], n: int, seed: int) -> list[SpanRead]:
    reads = synth.reads(synth.allele(inner), n, err=ERR, seed=seed)
    return categorize_reads(reads, ANCH, S).spanning


def _het(inner: list[str], seed_a: int, seed_b: int, n_b: int = N_PER_ALLELE) -> list[SpanRead]:
    return _spans(inner, N_PER_ALLELE, seed_a) + _spans(PLAIN, n_b, seed_b)


def _misgrouped(res: PhaseResult, first: list[SpanRead]) -> int:
    """Reads in a group whose majority origin differs from their own."""
    names = {m.name for m in first}
    wrong = 0
    for g in res.groups:
        own = sum(m.name in names for m in g)
        wrong += min(own, len(g) - own)
    return wrong


def _one_base_off_run() -> str:
    """X with the base right after a run changed to the run base (extends the run)."""
    x = synth.RD.repeats["X"]
    _start, e, b = next(
        r for r in _runs(x, S.phase_run_min_len) if r[1] + 1 < len(x) and x[r[1] + 1] != r[2]
    )
    return x[:e] + b + x[e + 1 :]


def _isolated(cons: str, start: int) -> int:
    """First position from start where no two neighbours in [p - 2, p + 2] are equal."""
    return next(
        p
        for p in range(start, len(cons))
        if all(cons[i] != cons[i + 1] for i in range(p - 2, p + 2))
    )


def _stranded(clean: str, edited: str, frac: float, seed: int) -> list[SpanRead]:
    """Homozygous sample whose "-" reads carry a systematic error at fraction frac."""
    n_edit = round(frac * N_PER_ALLELE)
    plus = synth.reads(clean, N_PER_ALLELE, err=ERR, seed=seed, strand_mix=False)
    minus = synth.reads(edited, n_edit, err=ERR, seed=seed + 1, strand_mix=False)
    minus += synth.reads(clean, N_PER_ALLELE - n_edit, err=ERR, seed=seed + 2, strand_mix=False)
    flipped = [ReadRecord(r.name, rc(r.seq), r.qual[::-1]) for r in minus]
    return categorize_reads(plus + flipped, ANCH, S).spanning


def test_equal_length_alleles_with_two_linked_differences_split() -> None:
    members = _het(TWO_UNITS, 1, 2)
    res = split_by_linked_sites(synth.allele(TWO_UNITS), members, S, random.Random(S.seed))
    assert res.basis == "linked_sites"
    assert sorted(len(g) for g in res.groups) == [N_PER_ALLELE, N_PER_ALLELE]
    assert len(res.sites) >= S.min_linked_sites and res.candidate is None
    assert res.unassigned == []


@pytest.mark.parametrize("inner", [TWO_UNITS, ONE_B], ids=["two_units", "one_b"])
@pytest.mark.parametrize(("seed_a", "seed_b"), PURITY_SEEDS)
def test_split_groups_are_pure_whatever_the_site_orientation(
    inner: list[str], seed_a: int, seed_b: int
) -> None:
    first = _spans(inner, N_PER_ALLELE, seed_a)
    members = first + _spans(PLAIN, N_PER_ALLELE, seed_b)
    res = split_by_linked_sites(synth.allele(inner), members, S, random.Random(S.seed))
    assert res.basis == "linked_sites" and _misgrouped(res, first) == 0
    # A read with conflicting event votes (a sequencing error at one of two sites) is
    # left unassigned rather than guessed; every member is either grouped or unassigned.
    assert sum(len(g) for g in res.groups) + len(res.unassigned) == len(members)
    names = {m.name for m in first}
    assert sorted(any(m.name in names for m in g) for g in res.groups) == [False, True]


def test_identical_alleles_do_not_split() -> None:
    members = _spans(PLAIN, 2 * N_PER_ALLELE, 3)
    res = split_by_linked_sites(synth.allele(PLAIN), members, S, random.Random(S.seed))
    assert res.basis == "none" and len(res.groups) == 1 and res.groups[0] == members
    assert res.unassigned == []


def test_single_base_unit_difference_is_unconfirmed() -> None:
    members = _het(ONE_Q, 1, 2)
    res = split_by_linked_sites(synth.allele(ONE_Q), members, S, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site" and len(res.groups) == 1
    assert res.candidate is not None and res.groups[0] == members


def test_multi_base_unit_difference_counts_as_linked_sites() -> None:
    members = _het(ONE_A, 1, 2)
    res = split_by_linked_sites(synth.allele(ONE_A), members, S, random.Random(S.seed))
    assert res.basis == "linked_sites"


def test_min_linked_sites_is_read_from_settings() -> None:
    # Linked events never outnumber linked sites, so one more than the site count can
    # never be confirmed.
    members = _het(TWO_UNITS, 1, 2)
    base = split_by_linked_sites(synth.allele(TWO_UNITS), members, S, random.Random(S.seed))
    strict = dataclasses.replace(S, min_linked_sites=len(base.sites) + 1)
    res = split_by_linked_sites(synth.allele(TWO_UNITS), members, strict, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site" and res.candidate is not None


def test_minority_group_below_het_min_group_is_not_split() -> None:
    # The largest minority m with m < het_min_group * (m + N) is linked (het_af_min and
    # the minor-read floor are loosened), but the group-size guard must refuse the split.
    # The minority is drawn evenly from both strands so it stays strand-consistent.
    members = _het(TWO_UNITS, 1, 2)
    minority = math.ceil(S.het_min_group * N_PER_ALLELE / (1 - S.het_min_group)) - 1
    first = members[:N_PER_ALLELE]
    plus = [m for m in first if m.strand == "+"][: (minority + 1) // 2]
    minus = [m for m in first if m.strand == "-"][: minority // 2]
    lopsided = plus + minus + members[N_PER_ALLELE:]
    loose = dataclasses.replace(
        S, het_af_min=minority / len(lopsided) / 2, phase_min_minor_reads=minority // 2
    )
    res = split_by_linked_sites(synth.allele(PLAIN), lopsided, loose, random.Random(S.seed))
    assert res.basis == "unconfirmed_group_size" and len(res.groups) == 1
    assert len(res.sites) >= S.min_linked_sites  # linked, but refused by the group guard
    top = max(res.sites, key=lambda site: site["af"] * (1 - site["af"]))
    assert res.candidate == top


def test_site_read_cap_subsamples_deterministically_through_rng() -> None:
    members = _het(TWO_UNITS, 1, 2)
    capped = dataclasses.replace(S, phase_max_site_reads=len(members) // 2)
    runs = [
        split_by_linked_sites(synth.allele(TWO_UNITS), members, capped, random.Random(S.seed))
        for _ in range(2)
    ]
    assert runs[0].sites == runs[1].sites and runs[0].basis == "linked_sites"
    # Grouping still covers every member, not only the sampled subset.
    assert sorted(len(g) for g in runs[0].groups) == [N_PER_ALLELE, N_PER_ALLELE]
    assert all(s["n"] == capped.phase_max_site_reads for s in runs[0].sites)


def test_minor_read_floor_is_read_from_settings() -> None:
    members = _het(TWO_UNITS, 1, 2)
    strict = dataclasses.replace(S, phase_min_minor_reads=N_PER_ALLELE + 1)
    res = split_by_linked_sites(synth.allele(TWO_UNITS), members, strict, random.Random(S.seed))
    assert res.basis == "none"


# Variant-draft seed pairs 3/4, 9/10 and 11/12 gave "none" while a run was its own
# background (before the leave-one-out fix); 1/2 did not.
@pytest.mark.parametrize(
    ("draft", "seed_a"),
    [("plain", 1), ("variant", 1), ("variant", 3), ("variant", 9), ("variant", 11)],
)
def test_snv_next_to_a_run_is_one_event(draft: str, seed_a: int) -> None:
    # Extending a run by one base changes both the run length and the adjacent column;
    # on the plain draft these touching features are one event, so the split is
    # unconfirmed. On the variant draft the run length is unique in the consensus, so
    # the leave-one-out run background must not hide it (it used to give "none").
    inner = ["X"] * 10 + [_one_base_off_run()] + ["X"] * 19
    members = _het(inner, seed_a, seed_a + 1)
    cons = synth.allele(PLAIN if draft == "plain" else inner)
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site" and len(res.groups) == 1
    kinds = {site["site"][0] for site in res.sites}
    assert kinds >= ({"run", "col"} if draft == "plain" else {"run"})


def _two_subs(cons: str) -> str:
    p1 = _isolated(cons, len(cons) // 2)
    p2 = _isolated(cons, p1 + synth.RD.repeat_length_bp)
    out = list(cons)
    for p in (p1, p2):
        out[p] = next(b for b in "ACGT" if b not in cons[p - 1 : p + 2])
    return "".join(out)


def _one_sub(cons: str) -> str:
    p = _isolated(cons, len(cons) // 2)
    return cons[:p] + next(b for b in "ACGT" if b not in cons[p - 1 : p + 2]) + cons[p + 1 :]


def _deletion(cons: str) -> str:
    p = _isolated(cons, len(cons) // 2)
    return cons[:p] + cons[p + 2 :]


@pytest.mark.parametrize(
    "edit", [_two_subs, _one_sub, _deletion], ids=["two_subs", "one_sub", "two_bp_deletion"]
)
@pytest.mark.parametrize("frac", STRAND_ERROR_FRACS)
def test_strand_specific_errors_do_not_split_a_homozygote(
    edit: Callable[[str], str], frac: float
) -> None:
    cons = synth.allele(PLAIN)
    members = _stranded(cons, edit(cons), frac, S.seed)
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    assert res.basis == "none" and len(res.groups) == 1


@pytest.mark.parametrize("inner", [ONE_B, TWO_UNITS, ONE_A], ids=["one_b", "two_units", "one_a"])
@pytest.mark.parametrize("ratio", ALLELE_RATIOS, ids=["3to1", "1to3"])
@pytest.mark.parametrize(("seed_a", "seed_b"), IMBALANCE_SEEDS)
def test_unbiased_imbalanced_heterozygote_still_splits(
    inner: list[str], ratio: tuple[int, int], seed_a: int, seed_b: int
) -> None:
    total = 2 * N_PER_ALLELE
    n_variant = total * ratio[0] // sum(ratio)
    first = _spans(inner, n_variant, seed_a)
    members = first + _spans(PLAIN, total - n_variant, seed_b)
    res = split_by_linked_sites(synth.allele(inner), members, S, random.Random(S.seed))
    assert res.basis == "linked_sites" and _misgrouped(res, first) == 0


@pytest.mark.parametrize("above", [False, True], ids=["noise_below_af_min", "noise_above"])
def test_run_with_no_same_base_peer_never_splits(above: bool) -> None:
    # A run of a base that forms no other run in the consensus has no background peer,
    # so its background falls back to 0 and only het_af_min applies. Run-length noise
    # there (one base lost) is a single event: at most an unconfirmed candidate.
    cons_plain = synth.allele(PLAIN)
    base = next(b for b in "ACGT" if all(r[2] != b for r in _runs(cons_plain, S.phase_run_min_len)))
    x = synth.RD.repeats["X"]
    pos = _isolated(x, S.phase_run_min_len)
    long_run = base * (S.phase_run_min_len + 1)
    unit, short = x[:pos] + long_run + x[pos:], x[:pos] + long_run[1:] + x[pos:]
    inner, noisy = ["X"] * 10 + [unit] + ["X"] * 19, ["X"] * 10 + [short] + ["X"] * 19
    total = 2 * N_PER_ALLELE
    frac = S.het_af_min * (S.phase_gap_af_factor if above else 1 / S.phase_gap_af_factor)
    n_noisy = round(frac * total)
    members = _spans(inner, total - n_noisy, S.seed) + _spans(noisy, n_noisy, S.seed + 1)
    cons = synth.allele(inner)
    runs = [r for r in _runs(cons, S.phase_run_min_len) if r[2] == base]
    assert len(runs) == 1  # the fallback-to-zero precondition
    res = split_by_linked_sites(cons, members, S, random.Random(S.seed))
    assert len(res.groups) == 1
    if above:
        assert res.basis == "unconfirmed_single_site"
        assert res.candidate is not None and res.candidate["site"] == ("run", runs[0][0])
    else:
        assert res.basis == "none"


def test_uninformative_and_tied_reads_are_unassigned_not_grouped() -> None:
    # ONE_B differs from PLAIN at two non-touching columns (two events). A read with a
    # third base at both is uninformative; one carrying B at the first and X at the
    # second ties between the groups. Both must be left for the caller to reassign.
    with_b, plain = synth.allele(ONE_B), synth.allele(PLAIN)
    first, second = (i for i, (x, y) in enumerate(zip(with_b, plain, strict=True)) if x != y)
    blank = list(with_b)
    for pos in (first, second):
        blank[pos] = next(b for b in "ACGT" if b not in (with_b[pos], plain[pos]))
    tie = with_b[:second] + plain[second] + with_b[second + 1 :]
    odd = [
        SpanRead("uninformative", "".join(blank), 0.0, "+", 0, "motif"),
        SpanRead("tie", tie, 0.0, "+", 0, "motif"),
    ]
    members = _het(ONE_B, 1, 2)
    res = split_by_linked_sites(with_b, [*members, *odd], S, random.Random(S.seed))
    assert res.basis == "linked_sites"
    assert odd[0] in res.unassigned and odd[1] in res.unassigned
    grouped = [m for g in res.groups for m in g]
    assert sorted(m.name for m in grouped + res.unassigned) == sorted(
        m.name for m in [*members, *odd]
    )


def test_pairwise_linkage_read_floor_is_read_from_settings() -> None:
    # No pair of sites can reach more informative reads than there are members, so
    # nothing links and the split stays unconfirmed.
    members = _het(TWO_UNITS, 1, 2)
    strict = dataclasses.replace(S, phase_min_pair_reads=len(members) + 1)
    res = split_by_linked_sites(synth.allele(TWO_UNITS), members, strict, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site" and res.candidate is not None
