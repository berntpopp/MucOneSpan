"""S2 smear significance test: exact conditional Poisson rate test and verdict bands."""

from __future__ import annotations

import math

import pytest

from muc_one_span.hybrid.smear import binomial_sf, smear_test, smear_verdict
from muc_one_span.settings import HybridSettings

S = HybridSettings()
UNIT = 60


@pytest.mark.parametrize(
    ("k", "n", "p", "expected"),
    [(0, 5, 0.3, 1.0), (6, 5, 0.3, 0.0), (1, 5, 0.0, 0.0), (1, 5, 1.0, 1.0)],
)
def test_binomial_sf_edges(k: int, n: int, p: float, expected: float) -> None:
    assert binomial_sf(k, n, p) == expected


def test_binomial_sf_matches_closed_form() -> None:
    n, p = 10, 0.3
    assert binomial_sf(n, n, p) == pytest.approx(p**n)
    assert binomial_sf(1, n, p) == pytest.approx(1 - (1 - p) ** n)
    pmf2 = math.comb(n, 2) * p**2 * (1 - p) ** (n - 2)
    assert binomial_sf(2, n, p) - binomial_sf(3, n, p) == pytest.approx(pmf2)


def test_isolated_cluster_with_no_background_is_highly_significant() -> None:
    lengths = [2340.0] * 5
    test = smear_test(2340.0, 50.0, lengths, (900.0, 4000.0), S, UNIT)
    assert test.core_reads == 5
    assert all(side.reads == 0 for side in test.sides)
    assert test.p_value < S.smear_test_alpha / S.smear_test_borderline_factor


def test_cluster_on_uniform_background_is_not_significant() -> None:
    lengths = [900.0 + 10.0 * i for i in range(311)]  # one read every 10 bp
    test = smear_test(2340.0, 50.0, lengths, (900.0, 4000.0), S, UNIT)
    assert test.p_value > S.smear_test_alpha * S.smear_test_borderline_factor


def test_side_clipped_to_zero_width_is_skipped() -> None:
    # The candidate's window reaches below the region floor: only the right side informs.
    test = smear_test(930.0, 50.0, [930.0] * 4, (900.0, 4000.0), S, UNIT)
    assert len(test.sides) == 1
    no_sides = smear_test(930.0, 50.0, [930.0] * 4, (900.0, 960.0), S, UNIT)
    assert (no_sides.sides, no_sides.p_value) == ((), 0.0)


def test_verdict_bands_and_correction() -> None:
    alpha, factor = S.smear_test_alpha, S.smear_test_borderline_factor
    assert smear_verdict(alpha / factor / 2, 1, S) == "significant"
    assert smear_verdict(alpha, 1, S) == "borderline"
    assert smear_verdict(alpha * factor, 1, S) == "smear"
    # Bonferroni over 4 candidates pushes a borderline p into 'smear'; "none" does not.
    assert smear_verdict(alpha, 4, S) == "smear"
    assert smear_verdict(alpha, 4, HybridSettings(smear_test_correction="none")) == "borderline"
