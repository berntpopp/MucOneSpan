"""Task 15f: run-length-aware stutter profiles (``hybrid.stutter``), exact inputs."""

from __future__ import annotations

import dataclasses
import math

import pytest

from muc_one_span.hybrid import stutter
from muc_one_span.settings import HybridSettings

S = HybridSettings()
SIZE = S.hp_max_run_len + 1


def _profile(length: int, shares: dict[int, float]) -> list[float]:
    """Observed-length profile of a run of ``length``: error -> share, rest exact."""
    p = [0.0] * SIZE
    for error, share in shares.items():
        p[length + error] = share
    p[length] = 1 - sum(shares.values())
    return p


def test_measured_length_is_used_directly() -> None:
    p7 = _profile(7, {-1: 0.08, 1: 0.05})
    assert stutter.length_profile({7: p7}, 7, S.hp_stutter_max_growth) == p7


def test_single_measured_length_is_shifted() -> None:
    p7 = _profile(7, {-1: 0.08, 1: 0.05})
    assert stutter.length_profile({7: p7}, 8, S.hp_stutter_max_growth) == stutter.shift(p7, 1)


def test_no_measured_length_gives_none() -> None:
    assert stutter.length_profile({}, 8, S.hp_stutter_max_growth) is None


def test_missing_length_is_extrapolated_per_error_value() -> None:
    """Each error share follows its own geometric trend from the two nearest lengths."""
    growth = S.hp_stutter_max_growth / 2
    base = {-1: 0.04, 1: 0.05}
    p6 = _profile(6, base)
    p7 = _profile(7, {-1: base[-1] * growth, 1: base[1]})
    p8 = stutter.length_profile({6: p6, 7: p7, 3: _profile(3, {})}, 8, S.hp_stutter_max_growth)
    assert p8 is not None
    assert math.isclose(sum(p8), 1.0)
    # Unnormalised: deletion x growth, insertion unchanged, exact reads x (p7/p6 exact).
    raw = {7: base[-1] * growth**2, 9: base[1], 8: p7[7] ** 2 / p6[6]}
    total = sum(raw.values())
    for k, v in raw.items():
        assert math.isclose(p8[k], v / total), (k, p8[k], v / total)


def test_extrapolated_growth_is_capped() -> None:
    cap = S.hp_stutter_max_growth
    p6 = _profile(6, {-1: 0.01})
    p7 = _profile(7, {-1: 0.01 * cap * cap})
    p8 = stutter.length_profile({6: p6, 7: p7}, 8, cap)
    assert p8 is not None
    exact = p7[7] * max(min(p7[7] / p6[6], cap), 1 / cap)
    total = p7[6] * cap + exact
    assert math.isclose(p8[7], p7[6] * cap / total)
    loose = stutter.length_profile({6: p6, 7: p7}, 8, cap * cap)
    assert loose is not None and loose[7] > p8[7]


def test_interpolation_between_measured_lengths() -> None:
    p6 = _profile(6, {-1: 0.02})
    p8 = _profile(8, {-1: 0.08})
    p7 = stutter.length_profile({6: p6, 8: p8}, 7, S.hp_stutter_max_growth)
    assert p7 is not None
    assert p6[5] < p7[6] < p8[7]


def test_growth_cap_is_validated() -> None:
    with pytest.raises(ValueError, match="hp_stutter_max_growth"):
        dataclasses.replace(S, hp_stutter_max_growth=0.99)
