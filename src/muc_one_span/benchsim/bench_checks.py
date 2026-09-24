"""Validation helpers and factor-level names shared by the bench config modules.

`bench_config` and `bench_sets` both validate with these helpers; they live
here so neither module imports the other's internals.
"""

from __future__ import annotations

import math

PCR_LEVEL_NAMES = ("calibrated", "strong", "none")  # semantics implemented in `profiles`
ERROR_LEVEL_NAMES = ("calibrated", "poor")
# Benchmark profiles (semantics in `design`, `profiles`, `generate`).
PROFILE_NAMES = ("ont_amplicon_r10", "ont_genomic_targeted", "hifi_amplicon")


def check_num(
    name: str, value: object, lo: float, hi: float | None = None, *, open_lo: bool = False
) -> None:
    """Raise unless ``value`` is a finite number in ``[lo, hi]`` (``(lo, hi]`` if open_lo)."""
    ok = False
    if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
        ok = (value > lo if open_lo else value >= lo) and (hi is None or value <= hi)
    if not ok:
        bound = f"{'(' if open_lo else '['}{lo}, {hi if hi is not None else 'inf'}]"
        raise ValueError(f"{name} must be a finite number in {bound}")


def check_int(name: str, value: object, lo: int) -> None:
    """Raise unless ``value`` is an integer (not a bool) >= ``lo``."""
    if type(value) is not int or value < lo:
        raise ValueError(f"{name} must be an integer >= {lo}")


def check_levels(name: str, values: object, lo: float, hi: float) -> None:
    """Raise unless ``values`` is a nonempty tuple of distinct numbers in ``[lo, hi]``."""
    if not isinstance(values, tuple) or not values:
        raise ValueError(f"{name} must be a nonempty list of distinct values")
    for value in values:
        check_num(name, value, lo, hi)
    if len(set(values)) != len(values):
        raise ValueError(f"{name} must be a nonempty list of distinct values")


def check_names(name: str, values: object, allowed: tuple[str, ...]) -> None:
    """Raise unless ``values`` is a nonempty tuple of distinct names from ``allowed``."""
    if not isinstance(values, tuple) or not values or not all(isinstance(v, str) for v in values):
        raise ValueError(f"{name} must be a nonempty list of distinct names")
    if len(set(values)) != len(values):
        raise ValueError(f"{name} must be a nonempty list of distinct names")
    if not set(values) <= set(allowed):
        raise ValueError(f"{name} must use only {allowed!r}")
