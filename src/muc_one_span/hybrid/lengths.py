"""S2: allele length peaks from spanning-read lengths with length-scaled windows."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

UNIT = 60
GATE_RELEVANT_REJECTIONS = frozenset({"support_below_threshold", "max_alleles"})


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


def window_bp(length_bp: float, settings: HybridSettings) -> float:
    """Assignment half-window: grows with allele length (ONT span noise ~0.25 bp/unit)."""
    return settings.peak_window_base_bp + settings.peak_window_per_unit_bp * length_bp / UNIT


def _density(lengths: list[float]) -> tuple[list[float], list[float]]:
    """Gaussian kernel density on a 2 bp grid; bandwidth grows with length."""
    lo, hi = min(lengths) - 100.0, max(lengths) + 100.0
    grid = [lo + 2.0 * i for i in range(int((hi - lo) / 2.0))]
    dens = [0.0] * len(grid)
    for value in lengths:
        bw = 8.0 + 0.004 * value
        first = max(0, int((value - 4 * bw - lo) / 2.0))
        last = min(len(grid), int((value + 4 * bw - lo) / 2.0) + 1)
        for i in range(first, last):
            dens[i] += math.exp(-0.5 * ((grid[i] - value) / bw) ** 2)
    return grid, dens


def _is_smear(c: float, top: float, lengths: list[float], settings: HybridSettings) -> bool:
    """Below the major peak by >1.5 units, a candidate is a short product (smear) unless it
    has at least ``min_peak_reads`` reads and stands out from its shoulders."""
    if c >= top - 1.5 * UNIT:
        return False
    w = window_bp(c, settings)
    inside = sum(abs(x - c) <= w for x in lengths)
    shoulders = sum(w < abs(x - c) <= 3 * w for x in lengths) / 2  # same width as inside
    return inside < settings.min_peak_reads or inside < settings.smear_min_prominence * max(
        shoulders, 1.0
    )


def _reason(
    c: float, top: float, support: int, n_kept: int, lengths: list[float], settings: HybridSettings
) -> str | None:
    """None when the candidate is accepted, else the rejection reason."""
    if _is_smear(c, top, lengths, settings):
        return "smear"
    if support <= settings.rejected_peak_noise_reads:
        return "noise"
    far = abs(c - top) >= 2 * UNIT
    frac = settings.far_peak_min_frac if far else settings.near_peak_min_frac
    top_support = sum(abs(x - top) <= window_bp(top, settings) for x in lengths)
    if support < settings.min_peak_reads or support < frac * top_support:
        return "support_below_threshold"
    return "max_alleles" if n_kept >= 2 else None


def fit_length_model(spans: list[SpanRead], settings: HybridSettings) -> LengthModel:
    """Pick up to two allele peaks; report rejected peaks, short products, unassigned reads."""
    if not spans:
        return LengthModel([], [], [], [], 0)
    lengths = [float(s.length) for s in spans]
    grid, dens = _density(lengths)
    maxima = [
        i for i in range(1, len(grid) - 1) if dens[i] >= dens[i - 1] and dens[i] > dens[i + 1]
    ]
    maxima.sort(key=lambda i: -dens[i])
    centers: list[float] = []
    for i in maxima:  # keep maxima at least ~0.7 unit apart (Δ1 alleles stay separable)
        if all(abs(grid[i] - c) >= 0.7 * UNIT for c in centers):
            centers.append(grid[i])
    support = {c: sum(abs(x - c) <= window_bp(c, settings) for x in lengths) for c in centers}
    top = max(centers, key=lambda c: support[c])
    kept, rejected = [top], []
    for c in sorted(centers, key=lambda c: -support[c]):
        if c == top:
            continue
        reason = _reason(c, top, support[c], len(kept), lengths, settings)
        if reason is None:
            kept.append(c)
            continue
        rejected.append(
            {
                "center_bp": round(c, 1),
                "units": round(c / UNIT),
                "support": support[c],
                "reason": reason,
            }
        )
    peaks = [LengthPeak(c, support[c]) for c in sorted(kept)]
    short_products: list[SpanRead] = []
    unassigned: list[SpanRead] = []
    shortest = peaks[0].center_bp
    for sp in spans:
        dist = [abs(sp.length - p.center_bp) for p in peaks]
        j = min(range(len(peaks)), key=lambda k: dist[k])
        if dist[j] <= window_bp(peaks[j].center_bp, settings):
            peaks[j].members.append(sp)
        elif sp.length < shortest - 1.5 * UNIT:
            short_products.append(sp)
        else:
            unassigned.append(sp)
    return LengthModel(peaks, rejected, short_products, unassigned, len(spans))
