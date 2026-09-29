"""S2 smear test: is a below-top length cluster a significant excess over local smear?

PCR/ONT smear (internally deleted products) spreads reads over a broad, uneven range of
lengths below the true allele. A chance cluster of smear reads grows with depth, so an
absolute read count cannot separate it from a real minor allele. Instead, each
candidate's read count in a narrow *core* window is compared with the local smear
background on each side of its assignment window, using the exact conditional test for
two Poisson rates: given ``k`` core reads and ``m`` background reads, ``k`` is
Binomial(k + m, rho / (1 + rho)) under "same density", where ``rho`` is the core width
over the background width. Because the background count enters the test, the
uncertainty of the background estimate is accounted for (a plug-in Poisson test with an
estimated rate is anti-conservative).

The candidate must exceed the background on *both* sides (the larger of the two p
values is used), so a cluster at the edge of a steep smear density step is not called
significant against the empty side. Sides clipped to zero width by the region carry no
information and are skipped; with no informative side the p value is 0 (nothing
explains the cluster as smear).
"""

from __future__ import annotations

import math
from typing import NamedTuple

from muc_one_span.settings import HybridSettings


class SideBackground(NamedTuple):
    """Reads counted on one side of the candidate and the width they were counted over."""

    reads: int
    width_bp: float


class SmearTest(NamedTuple):
    """Core count, per-side backgrounds, and the uncorrected one-sided p value."""

    core_reads: int
    core_half_bp: float
    sides: tuple[SideBackground, ...]
    p_value: float


def binomial_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p), summed in log space (pure Python)."""
    if k <= 0:
        return 1.0
    if k > n or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    log_p, log_q, log_n = math.log(p), math.log1p(-p), math.lgamma(n + 1)
    total = sum(
        math.exp(log_n - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * log_p + (n - i) * log_q)
        for i in range(k, n + 1)
    )
    return min(1.0, total)


def _side(
    distances: list[float], limit_bp: float, settings: HybridSettings, unit_bp: int
) -> SideBackground | None:
    """Background on one side: at least the configured flank, widened until it holds
    ``smear_background_min_reads`` reads, never past the region edge (``limit_bp``)."""
    if limit_bp <= 0:
        return None
    distances = sorted(distances)
    need = settings.smear_background_min_reads
    flank = settings.smear_background_flank_units * unit_bp
    # Too few reads anywhere on this side: the whole side is the background.
    reach = max(flank, distances[need - 1]) if len(distances) >= need else limit_bp
    reach = min(reach, limit_bp)
    return SideBackground(sum(1 for d in distances if d <= reach), reach)


def smear_test(
    center_bp: float,
    window_half_bp: float,
    lengths: list[float],
    region: tuple[float, float],
    settings: HybridSettings,
    unit_bp: int,
) -> SmearTest:
    """One-sided test of the candidate's core count against its local smear background.

    ``region`` is the below-top length range ``(lo, hi)`` where smear is assessed; the
    background sides start at the candidate's assignment-window edge
    (``window_half_bp``) so the candidate's own reads never inflate its background.
    """
    lo, hi = region
    core_half = settings.smear_test_window_frac * window_half_bp
    core = sum(abs(x - center_bp) <= core_half for x in lengths)
    left_edge, right_edge = center_bp - window_half_bp, center_bp + window_half_bp
    left = [left_edge - x for x in lengths if lo <= x < left_edge]
    right = [x - right_edge for x in lengths if right_edge < x <= hi]
    sides = tuple(
        side
        for side in (
            _side(left, left_edge - lo, settings, unit_bp),
            _side(right, hi - right_edge, settings, unit_bp),
        )
        if side is not None
    )
    p_value = 0.0
    for side in sides:
        rho = 2 * core_half / side.width_bp
        p_value = max(p_value, binomial_sf(core, core + side.reads, rho / (1 + rho)))
    return SmearTest(core, core_half, sides, p_value)


def smear_verdict(p_value: float, n_tests: int, settings: HybridSettings) -> str:
    """'significant', 'borderline' or 'smear' for one candidate's corrected p value.

    ``smear_test_correction="bonferroni"`` multiplies by the number of below-top
    candidates tested. The borderline band is ``[alpha / factor, alpha * factor)``.
    """
    family = n_tests if settings.smear_test_correction == "bonferroni" else 1
    adjusted = min(1.0, p_value * max(family, 1))
    alpha, factor = settings.smear_test_alpha, settings.smear_test_borderline_factor
    if adjusted < alpha / factor:
        return "significant"
    if adjusted < alpha * factor:
        return "borderline"
    return "smear"
