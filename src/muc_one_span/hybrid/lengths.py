"""S2: allele length peaks from spanning-read lengths with length-scaled windows.

``unit_bp`` (the MUC1 repeat-unit length) is never hardcoded here: callers derive it
from the loaded repeat dictionary (``Anchors.unit_bp``, see ``hybrid.spans``) and pass
it in explicitly, so this module has no biological or tuning constant of its own.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, NamedTuple

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


class _Background(NamedTuple):
    """Local window occupancy plus the wider below-top region around a candidate."""

    inside: int
    window_bp: float
    outside: int
    region_width: float


def _region_background(
    c: float, top: float, lengths: list[float], settings: HybridSettings, unit_bp: int
) -> _Background:
    """Read counts inside the candidate's own window versus the rest of the below-top
    region. ``outside == 0`` means there is no smear/background context at all around
    this candidate -- it is an isolated cluster, not debris in a smear field (I2)."""
    w = window_bp(c, settings, unit_bp)
    inside = _count_within(lengths, c, w)
    region_lo = min(lengths)
    region_hi = top - settings.smear_short_product_units * unit_bp
    region_width = max(region_hi - region_lo, 1.0)
    total_in_region = sum(1 for x in lengths if region_lo <= x <= region_hi)
    outside = max(total_in_region - inside, 0)
    return _Background(inside, w, outside, region_width)


def _is_smear(lengths: list[float], c: float, settings: HybridSettings, bg: _Background) -> bool:
    """A candidate with background around it is smear unless it has at least
    ``min_peak_reads`` reads and stands out from its local shoulders."""
    mult = settings.smear_shoulder_width_mult
    shoulders = sum(bg.window_bp < abs(x - c) <= mult * bg.window_bp for x in lengths) / (mult - 1)
    return bg.inside < settings.min_peak_reads or bg.inside < settings.smear_min_prominence * max(
        shoulders, settings.smear_shoulder_floor
    )


def _is_smear_ambiguous(settings: HybridSettings, bg: _Background) -> bool:
    """A candidate that clears the local shoulder check can still be a chance density bump
    in smear debris: this compares it against the *whole* below-top background rate rather
    than just its narrow shoulders, catching the cases the local check misses."""
    expected = (bg.outside / bg.region_width) * (2 * bg.window_bp)
    ratio = bg.inside / max(expected, settings.smear_background_floor)
    return ratio < settings.smear_background_ratio_min


def _reason(
    c: float,
    top: float,
    support: int,
    n_kept: int,
    lengths: list[float],
    settings: HybridSettings,
    unit_bp: int,
) -> str | None:
    """None when the candidate is accepted, else the rejection reason."""
    if c < top - settings.smear_short_product_units * unit_bp:
        bg = _region_background(c, top, lengths, settings, unit_bp)
        if bg.outside > 0:
            if _is_smear(lengths, c, settings, bg):
                return "smear"
            if _is_smear_ambiguous(settings, bg):
                return "smear_ambiguous"
    if support <= settings.rejected_peak_noise_reads:
        return "noise"
    far = abs(c - top) >= settings.peak_far_near_boundary_units * unit_bp
    frac = settings.far_peak_min_frac if far else settings.near_peak_min_frac
    top_support = _count_within(lengths, top, window_bp(top, settings, unit_bp))
    if support < settings.min_peak_reads or support < frac * top_support:
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
    kept, rejected = [top], []
    for c in sorted(centers, key=lambda c: -support[c]):
        if c == top:
            continue
        reason = _reason(c, top, support[c], len(kept), lengths, settings, unit_bp)
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
