"""S2 length peaks: scaled windows, rejected peaks, smear and close alleles."""

from __future__ import annotations

from functools import cache

from muc_one_span.hybrid.lengths import GATE_RELEVANT_REJECTIONS, LengthModel, fit_length_model
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S)
UNIT = ANCH.unit_bp


def spans(inner_units: int, n: int, seed: int, **kw: float) -> list[SpanRead]:
    seq = synth.allele(["X"] * inner_units)
    return categorize_reads(synth.reads(seq, n, err=0.02, seed=seed, **kw), ANCH, S).spanning


def _fit(reads: list[SpanRead]) -> LengthModel:
    return fit_length_model(reads, S, UNIT)


def test_two_distant_alleles_with_minor_long_allele() -> None:
    model = _fit(spans(30, 200, 1) + spans(70, 12, 2))
    assert sorted(round(p.center_bp / UNIT) for p in model.peaks) == [39, 79]
    assert model.gate_relevant_rejections == []


def test_smear_is_short_product_not_allele() -> None:
    model = _fit(spans(60, 120, 3, smear_frac=0.45))
    assert len(model.peaks) == 1
    assert model.short_product_fraction > 0.3
    assert model.gate_relevant_rejections == []


def test_real_short_allele_amid_smear_is_kept() -> None:
    model = _fit(spans(30, 40, 8) + spans(60, 120, 9, smear_frac=0.3))
    assert sorted(round(p.center_bp / UNIT) for p in model.peaks) == [39, 69]


def test_close_alleles_are_not_silently_dropped() -> None:
    # Alleles one unit apart must yield two peaks or a large unassigned fraction.
    model = _fit(spans(40, 150, 4) + spans(41, 150, 5))
    assert len(model.peaks) == 2 or model.unassigned_fraction > 0.2


def test_rejected_minor_peak_is_reported() -> None:
    model = _fit(spans(30, 200, 6) + spans(80, 5, 7))
    assert len(model.peaks) == 1
    assert [(r["support"], r["reason"]) for r in model.gate_relevant_rejections] == [
        (5, "support_below_threshold")
    ]


def test_third_real_peak_is_max_alleles() -> None:
    model = _fit(spans(30, 60, 10) + spans(50, 60, 11) + spans(70, 60, 12))
    assert len(model.peaks) == 2
    assert [r["reason"] for r in model.gate_relevant_rejections] == ["max_alleles"]


def test_single_outlier_read_is_noise_not_gate_relevant() -> None:
    model = _fit(spans(30, 100, 13) + spans(90, 1, 14))
    assert [r["reason"] for r in model.rejected] == ["noise"]
    assert model.gate_relevant_rejections == []


def test_empty_spans_returns_an_empty_model() -> None:
    model = _fit([])
    assert (model.peaks, model.rejected, model.short_products, model.unassigned) == ([], [], [], [])
    assert model.total == 0
    assert model.unassigned_fraction == 0.0
    assert model.short_product_fraction == 0.0


def test_isolated_short_allele_with_no_smear_background_is_support_below_threshold() -> None:
    # I2 (fix round 1): a below-top cluster with no reads anywhere else in the region (no
    # smear_frac, so no smear reads at all) must never be called "smear" or "smear_ambiguous";
    # a real-but-weak minor allele is a gate-relevant support call, not a silently discarded
    # smear artefact.
    model = _fit(spans(80, 200, 6) + spans(30, 5, 7))
    assert len(model.peaks) == 1
    assert [(r["support"], r["reason"]) for r in model.rejected] == [(5, "support_below_threshold")]
    assert model.gate_relevant_rejections == model.rejected


# Fixture geometry: every "minor" below is inner_units=30, bracketed by the same PRE/POST
# flank repeats as every other allele this file builds (synth.allele/synth.PRE/POST), so
# its KDE peak always falls at this many total repeat units.
MINOR_INNER_UNITS = 30
MINOR_UNITS = len(synth.PRE) + MINOR_INNER_UNITS + len(synth.POST)
MAJOR_INNER_UNITS = 60

# Fix round 4 acceptance grid (controller ruling: statistical significance of the local
# excess). Every criterion is asserted in EVERY cell; there is no exception list.
DS = (60, 120, 200, 300, 1000)
FRACS = (0.3, 0.45, 0.54)
GRID_SEEDS = 10  # seeds 0..9 per cell; runtime is reported in task-5-report.md
MINOR_SEED_OFFSET = 1000
F2_SEEDS = 8
MAX_HOMOZYGOUS_GATED_RATE = 0.10  # criterion (b)
MIN_MINOR_ACCEPT_RATE = 0.90  # criterion (c)
MIN_WEAK_MINOR_FLAGGED_RATE = 0.90  # criterion (c')
MINOR_SUPPORT_FRAC_20 = 0.2  # (c): genuine minor with 20% of the major's reads
MINOR_SUPPORT_FRAC_10 = 0.1  # (c'): weak minor with 10% of the major's reads
WEAK_MINOR_MIN_DEPTH = 120  # (c') applies at D >= 120


@cache
def _cached_spans(inner_units: int, n: int, seed: int, smear_frac: float) -> tuple[SpanRead, ...]:
    """Grid draws are shared by the (a)/(b)/(c)/(c') tests; generate each only once."""
    return tuple(spans(inner_units, n, seed, smear_frac=smear_frac))


def _major(depth: int, frac: float, seed: int) -> list[SpanRead]:
    return list(_cached_spans(MAJOR_INNER_UNITS, depth, seed, frac))


def _with_minor(depth: int, frac: float, seed: int, minor_frac: float) -> LengthModel:
    minor_n = round(minor_frac * depth)
    minor = list(_cached_spans(MINOR_INNER_UNITS, minor_n, seed + MINOR_SEED_OFFSET, 0.0))
    return _fit(_major(depth, frac, seed) + minor)


def _minor_accepted(model: LengthModel) -> bool:
    return MINOR_UNITS in {round(p.center_bp / UNIT) for p in model.peaks}


def _minor_flagged(model: LengthModel) -> bool:
    return any(
        r["units"] == MINOR_UNITS and r["reason"] in GATE_RELEVANT_REJECTIONS
        for r in model.rejected
    )


def test_homozygous_smear_never_silently_accepts_and_is_rarely_gate_relevant() -> None:
    """(a) a single allele plus smear never yields an unflagged second peak (zero
    tolerance); (b) it is gate-relevant (-> INCONCLUSIVE) in at most 10% of seeds, in
    every (D, smear_frac) cell up to D=1000."""
    false_accepts = []
    gated_rates = {}
    for depth in DS:
        for frac in FRACS:
            gated = 0
            for seed in range(GRID_SEEDS):
                model = _fit(_major(depth, frac, seed))
                if len(model.peaks) == 2 and model.gate_relevant_rejections == []:
                    false_accepts.append((depth, frac, seed))
                gated += bool(model.gate_relevant_rejections)
            gated_rates[(depth, frac)] = gated / GRID_SEEDS
    assert false_accepts == []
    over = {cell: rate for cell, rate in gated_rates.items() if rate > MAX_HOMOZYGOUS_GATED_RATE}
    assert over == {}


def test_genuine_twenty_percent_minor_is_accepted() -> None:
    """(c) a below-top minor with 20% of the major's reads is accepted in >= 90% of
    seeds, in every cell."""
    rates = {}
    for depth in DS:
        for frac in FRACS:
            hits = sum(
                _minor_accepted(_with_minor(depth, frac, seed, MINOR_SUPPORT_FRAC_20))
                for seed in range(GRID_SEEDS)
            )
            rates[(depth, frac)] = hits / GRID_SEEDS
    low = {cell: rate for cell, rate in rates.items() if rate < MIN_MINOR_ACCEPT_RATE}
    assert low == {}


def test_weak_ten_percent_minor_is_accepted_or_gate_relevant() -> None:
    """(c') a minor with 10% of the major's reads at D >= 120 is accepted or
    gate-relevant (never silent 'smear') in >= 90% of seeds, in every cell."""
    rates = {}
    for depth in (d for d in DS if d >= WEAK_MINOR_MIN_DEPTH):
        for frac in FRACS:
            hits = 0
            for seed in range(GRID_SEEDS):
                model = _with_minor(depth, frac, seed, MINOR_SUPPORT_FRAC_10)
                hits += _minor_accepted(model) or _minor_flagged(model)
            rates[(depth, frac)] = hits / GRID_SEEDS
    low = {cell: rate for cell, rate in rates.items() if rate < MIN_WEAK_MINOR_FLAGGED_RATE}
    assert low == {}


def test_isolated_minor_with_light_smear_is_support_below_threshold() -> None:
    # (e), F2: an isolated real 5-read minor stays support_below_threshold (gate-relevant),
    # never silently "smear", across light smear fractions (1-10% of the major smeared).
    for smear_frac in (0.01, 0.02, 0.05, 0.1):
        for seed in range(F2_SEEDS):
            model = _fit(spans(80, 200, 6 + seed, smear_frac=smear_frac) + spans(30, 5, 7 + seed))
            reasons = [r["reason"] for r in model.rejected if r["units"] == MINOR_UNITS]
            assert reasons == ["support_below_threshold"], (smear_frac, seed, model.rejected)


def test_clean_isolated_minor_is_accepted_not_ambiguous() -> None:
    """R1 (fix round 3): a below-top candidate with no smear background at all, that
    already clears the allele threshold, is accepted -- being isolated is not itself
    ambiguous."""
    for depth, minor_n in ((200, 50), (200, 20), (120, 30), (60, 15)):
        for seed in range(10):
            model = _fit(spans(80, depth, seed) + spans(MINOR_INNER_UNITS, minor_n, seed + 7777))
            assert _minor_accepted(model), (depth, minor_n, seed, model.rejected)


def test_borderline_smear_candidate_is_smear_ambiguous() -> None:
    # Widening the borderline band to cover every non-trivial p value routes the smear
    # tail's tested candidates to smear_ambiguous (gate-relevant) instead of silent smear.
    wide = HybridSettings(smear_test_borderline_factor=1e6)
    model = fit_length_model(spans(60, 120, 3, smear_frac=0.45), wide, UNIT)
    assert "smear_ambiguous" in [r["reason"] for r in model.gate_relevant_rejections]
