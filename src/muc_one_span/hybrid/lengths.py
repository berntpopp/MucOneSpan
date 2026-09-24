"""S2: allele length peaks from spanning-read lengths with length-scaled windows.

``unit_bp`` (the MUC1 repeat-unit length) is never hardcoded here: callers derive it
from the loaded repeat dictionary (``Anchors.unit_bp``, see ``hybrid.spans``) and pass
it in explicitly, so this module has no biological or tuning constant of its own.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from muc_one_span.hybrid.smear import smear_test, smear_verdict
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

GATE_RELEVANT_REJECTIONS = frozenset({"support_below_threshold", "max_alleles", "smear_ambiguous"})


@dataclass
class LengthPeak:
    """One accepted allele length peak and its spanning members."""

    center_bp: float
    support: int
    members: list[SpanRead] = field(default_factory=list)


@dataclass
class LengthModel:
    """Accepted peaks plus everything that was not used, so nothing is dropped silently."""

    peaks: list[LengthPeak]
    rejected: list[dict[str, Any]]
    short_products: list[SpanRead]
    unassigned: list[SpanRead]
    total: int

    @property
    def unassigned_fraction(self) -> float:
        return len(self.unassigned) / self.total if self.total else 0.0

    @property
    def short_product_fraction(self) -> float:
        return len(self.short_products) / self.total if self.total else 0.0

    @property
    def gate_relevant_rejections(self) -> list[dict[str, Any]]:
        """Rejected candidates that may be a real allele (not smear or single-read noise)."""
        return [r for r in self.rejected if r["reason"] in GATE_RELEVANT_REJECTIONS]


def window_bp(length_bp: float, settings: HybridSettings, unit_bp: int) -> float:
    """Assignment half-window: grows with allele length (ONT span noise ~0.25 bp/unit)."""
    return settings.peak_window_base_bp + settings.peak_window_per_unit_bp * length_bp / unit_bp


def _density(lengths: list[float], settings: HybridSettings) -> tuple[list[float], list[float]]:
    """Gaussian kernel density on a settings-scaled grid; bandwidth grows with length."""
    margin = settings.kde_grid_margin_bp
    step = settings.kde_grid_step_bp
    lo, hi = min(lengths) - margin, max(lengths) + margin
    grid = [lo + step * i for i in range(int((hi - lo) / step))]
    dens = [0.0] * len(grid)
    for value in lengths:
        bw = settings.kde_bandwidth_base_bp + settings.kde_bandwidth_per_bp * value
        reach = settings.kde_kernel_truncation_bw * bw
        first = max(0, int((value - reach - lo) / step))
        last = min(len(grid), int((value + reach - lo) / step) + 1)
        for i in range(first, last):
            # -0.5 * z**2 is the fixed Gaussian kernel exponent (defines the kernel
            # shape, not a tunable), so it stays a literal per the normal-density formula.
            dens[i] += math.exp(-0.5 * ((grid[i] - value) / bw) ** 2)
    return grid, dens


def _count_within(lengths: list[float], c: float, w: float) -> int:
    return sum(abs(x - c) <= w for x in lengths)


def _reason(
    c: float,
    top: float,
    support: int,
    n_kept: int,
    lengths: list[float],
    settings: HybridSettings,
    unit_bp: int,
    smear_ctx: tuple[tuple[float, float], int],
) -> str | None:
    """None when the candidate is accepted, else the rejection reason.

    Below-top candidates (more than ``smear_short_product_units`` below the top peak)
    first face the smear significance test (``hybrid.smear``, C4.2 fix round 4):
    not significant -> silent 'smear', whatever the absolute read count (smear debris
    grows with depth); in the configured borderline band around alpha ->
    'smear_ambiguous' (gate-relevant). Significant candidates, and every candidate that
    is not below-top, then need ``support >= max(min_peak_reads, frac * n_total)`` --
    total depth, not the top peak's own support (spec S2: "n_min, f_far/near * N") --
    else 'support_below_threshold' (gate-relevant). ``smear_ctx`` is the below-top
    region ``(lo, hi)`` and the number of candidates tested there (the correction family).
    """
    if support <= settings.rejected_peak_noise_reads:
        return "noise"
    region, n_tests = smear_ctx
    if c < region[1]:
        test = smear_test(c, window_bp(c, settings, unit_bp), lengths, region, settings, unit_bp)
        verdict = smear_verdict(test.p_value, n_tests, settings)
        if verdict == "smear":
            return "smear"
        if verdict == "borderline":
            return "smear_ambiguous"
    far = abs(c - top) >= settings.peak_far_near_boundary_units * unit_bp
    frac = settings.far_peak_min_frac if far else settings.near_peak_min_frac
    if support < max(settings.min_peak_reads, frac * len(lengths)):
        return "support_below_threshold"
    return "max_alleles" if n_kept >= 2 else None


def fit_length_model(spans: list[SpanRead], settings: HybridSettings, unit_bp: int) -> LengthModel:
    """Pick up to two allele peaks; report rejected peaks, short products, unassigned reads."""
    if not spans:
        return LengthModel([], [], [], [], 0)
    lengths = [float(s.length) for s in spans]
    grid, dens = _density(lengths, settings)
    maxima = [
        i for i in range(1, len(grid) - 1) if dens[i] >= dens[i - 1] and dens[i] > dens[i + 1]
    ]
    maxima.sort(key=lambda i: -dens[i])
    centers: list[float] = []
    min_sep = settings.peak_min_separation_units * unit_bp
    for i in maxima:  # keep maxima at least ~peak_min_separation_units apart (stay separable)
        if all(abs(grid[i] - c) >= min_sep for c in centers):
            centers.append(grid[i])
    support = {c: _count_within(lengths, c, window_bp(c, settings, unit_bp)) for c in centers}
    top = max(centers, key=lambda c: support[c])
    # Below-top smear region: from the shortest observable span to just below the top peak.
    region = (
        min(min(lengths), settings.min_span_units * unit_bp),
        top - settings.smear_short_product_units * unit_bp,
    )
    n_tests = sum(
        1
        for c in centers
        if c != top and c < region[1] and support[c] > settings.rejected_peak_noise_reads
    )
    kept, rejected = [top], []
    for c in sorted(centers, key=lambda c: -support[c]):
        if c == top:
            continue
        reason = _reason(
            c, top, support[c], len(kept), lengths, settings, unit_bp, (region, n_tests)
        )
        if reason is None:
            kept.append(c)
            continue
        rejected.append(
            {
                "center_bp": round(c, 1),
                "units": round(c / unit_bp),
                "support": support[c],
                "reason": reason,
            }
        )
    peaks = [LengthPeak(c, support[c]) for c in sorted(kept)]
    short_products: list[SpanRead] = []
    unassigned: list[SpanRead] = []
    shortest = peaks[0].center_bp
    short_cutoff = settings.smear_short_product_units * unit_bp
    for sp in spans:
        dist = [abs(sp.length - p.center_bp) for p in peaks]
        j = min(range(len(peaks)), key=lambda k: dist[k])
        if dist[j] <= window_bp(peaks[j].center_bp, settings, unit_bp):
            peaks[j].members.append(sp)
        elif sp.length < shortest - short_cutoff:
            short_products.append(sp)
        else:
            unassigned.append(sp)
    return LengthModel(peaks, rejected, short_products, unassigned, len(spans))
