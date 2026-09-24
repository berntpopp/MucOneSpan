"""S4 linked-site phase split for equal and close allele lengths."""

from __future__ import annotations

import dataclasses
import math
import random

import pytest

from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")

# Every phase tunable comes from S (HybridSettings defaults); nothing is duplicated here
# as a bare literal (owner directive: fully config-driven, no magic numbers).
S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S)
N_PER_ALLELE = 40
ERR = 0.02
PLAIN = ["X"] * 30
# A and B differ from X at 4 and 2 bases; Q differs from X at a single base.
TWO_UNITS = ["X"] * 10 + ["A"] + ["X"] * 10 + ["B"] + ["X"] * 8
ONE_A = ["X"] * 10 + ["A"] + ["X"] * 19
ONE_Q = ["X"] * 10 + ["Q"] + ["X"] * 19


def _spans(inner: list[str], n: int, seed: int) -> list[SpanRead]:
    reads = synth.reads(synth.allele(inner), n, err=ERR, seed=seed)
    return categorize_reads(reads, ANCH, S).spanning


def _het(inner: list[str], seed_a: int, seed_b: int, n_b: int = N_PER_ALLELE) -> list[SpanRead]:
    return _spans(inner, N_PER_ALLELE, seed_a) + _spans(PLAIN, n_b, seed_b)


def test_equal_length_alleles_with_two_linked_differences_split() -> None:
    members = _het(TWO_UNITS, 1, 2)
    res = split_by_linked_sites(synth.allele(TWO_UNITS), members, S, random.Random(S.seed))
    assert res.basis == "linked_sites"
    assert sorted(len(g) for g in res.groups) == [N_PER_ALLELE, N_PER_ALLELE]
    assert len(res.sites) >= S.min_linked_sites and res.candidate is None


@pytest.mark.parametrize(("seed_a", "seed_b"), [(5, 6), (7, 8), (9, 10)])
def test_split_groups_match_the_allele_of_origin(seed_a: int, seed_b: int) -> None:
    first = _spans(TWO_UNITS, N_PER_ALLELE, seed_a)
    members = first + _spans(PLAIN, N_PER_ALLELE, seed_b)
    res = split_by_linked_sites(synth.allele(TWO_UNITS), members, S, random.Random(S.seed))
    assert res.basis == "linked_sites"
    names = {m.name for m in first}
    assert sorted({m.name in names for m in g}.pop() for g in res.groups) == [False, True]
    assert all(len({m.name in names for m in g}) == 1 for g in res.groups)


def test_identical_alleles_do_not_split() -> None:
    members = _spans(PLAIN, 2 * N_PER_ALLELE, 3)
    res = split_by_linked_sites(synth.allele(PLAIN), members, S, random.Random(S.seed))
    assert res.basis == "none" and len(res.groups) == 1 and res.groups[0] == members


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
    # B differs from X at 2 bases and A at 4, so TWO_UNITS has 6 linked sites; one more
    # than that can never be confirmed.
    members = _het(TWO_UNITS, 1, 2)
    base = split_by_linked_sites(synth.allele(TWO_UNITS), members, S, random.Random(S.seed))
    strict = dataclasses.replace(S, min_linked_sites=len(base.sites) + 1)
    res = split_by_linked_sites(synth.allele(TWO_UNITS), members, strict, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site" and res.candidate is not None


def test_minority_group_below_het_min_group_is_not_split() -> None:
    # The largest minority m with m < het_min_group * (m + N) is linked (het_af_min and
    # the minor-read floor are loosened), but the group-size guard must refuse the split.
    members = _het(TWO_UNITS, 1, 2)
    minority = math.ceil(S.het_min_group * N_PER_ALLELE / (1 - S.het_min_group)) - 1
    lopsided = members[:minority] + members[N_PER_ALLELE:]
    loose = dataclasses.replace(
        S, het_af_min=minority / len(lopsided) / 2, phase_min_minor_reads=minority
    )
    res = split_by_linked_sites(synth.allele(PLAIN), lopsided, loose, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site" and len(res.groups) == 1
    assert len(res.sites) >= S.min_linked_sites  # linked, but refused by the group guard


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
