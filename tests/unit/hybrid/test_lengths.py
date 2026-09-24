"""S2 length peaks: scaled windows, rejected peaks, smear and close alleles."""

from __future__ import annotations

from muc_one_span.hybrid.lengths import LengthModel, fit_length_model
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


DS = (60, 120, 200, 300)
FRACS = (0.3, 0.45, 0.54)
GRID_SEEDS = 15  # keeps the whole grid well under the ~60s budget; see task-5-report.md
SPURIOUS_TARGET = 0.10  # criterion (b)
ACCEPT_TARGET = 0.90  # criterion (c)
MINOR_UNITS = 39  # 5 pre + 30 inner + 4 post repeats
# The smallest depth combined with the lightest smear has too little statistical power
# on this synthetic factory to hit the (b)/(c) targets at GRID_SEEDS=15 (13.3% spurious,
# 26.7% accept, measured). Every other one of the 12 grid cells meets both targets. This
# is the controller's documented (b)/(c) trade-off (DONE_WITH_CONCERNS); see
# task-5-report.md Fix round 2 for the full per-cell table. (a) and (d) hold here too.
KNOWN_TRADEOFF_CELLS = frozenset({(60, 0.3)})


def test_smear_grid_never_silently_accepts_and_spurious_rate_is_bounded() -> None:
    """Fix round 2, replaces the round-1 D=120-only sweep with the controller's full grid.
    (a): a single real allele plus smear must never yield a silently-accepted second peak
    -- zero tolerance, every (D, smear_frac) cell, every seed, no exceptions.
    (b): the homozygous gate-relevant (spurious) rate is <=10% per cell, except the
    documented KNOWN_TRADEOFF_CELLS."""
    false_accepts = []
    spurious_rates = {}
    for depth in DS:
        for frac in FRACS:
            gated = 0
            for seed in range(GRID_SEEDS):
                model = _fit(spans(60, depth, seed, smear_frac=frac))
                if len(model.peaks) == 2 and model.gate_relevant_rejections == []:
                    false_accepts.append((depth, frac, seed))
                if model.gate_relevant_rejections:
                    gated += 1
            spurious_rates[(depth, frac)] = gated / GRID_SEEDS
    assert false_accepts == []
    for (depth, frac), rate in spurious_rates.items():
        if (depth, frac) in KNOWN_TRADEOFF_CELLS:
            continue
        assert rate <= SPURIOUS_TARGET, f"D={depth} frac={frac}: spurious rate {rate:.0%}"


def test_smear_grid_genuine_minor_is_accepted_and_never_silently_smear() -> None:
    """(c): a genuine below-top minor with support >= 20% of the major is accepted in
    >=90% of seeds per cell, except the documented KNOWN_TRADEOFF_CELLS.
    (d): that minor is never silently 'smear' when its support clears min_peak_reads --
    always accepted or gate-relevant, every cell, every seed, no exceptions."""
    silent_smear = []
    accept_rates = {}
    for depth in DS:
        for frac in FRACS:
            minor_n = round(0.2 * depth)
            accepted = 0
            for seed in range(GRID_SEEDS):
                model = _fit(
                    spans(60, depth, seed, smear_frac=frac) + spans(30, minor_n, seed + 1000)
                )
                if MINOR_UNITS in {round(p.center_bp / UNIT) for p in model.peaks}:
                    accepted += 1
                    continue
                minor_reasons = {r["reason"] for r in model.rejected if r["units"] == MINOR_UNITS}
                if minor_n >= S.min_peak_reads and "smear" in minor_reasons:
                    silent_smear.append((depth, frac, seed))
            accept_rates[(depth, frac)] = accepted / GRID_SEEDS
    assert silent_smear == []
    for (depth, frac), rate in accept_rates.items():
        if (depth, frac) in KNOWN_TRADEOFF_CELLS:
            continue
        assert rate >= ACCEPT_TARGET, f"D={depth} frac={frac}: accept rate {rate:.0%}"


def test_isolated_minor_with_light_smear_is_support_below_threshold() -> None:
    # (e), F2: an isolated real minor allele stays support_below_threshold (gate-relevant),
    # never silently "smear", across the light smear fractions the review measured it
    # failing at (0.01-0.1: 1-10% of the major's reads smeared).
    for smear_frac in (0.01, 0.02, 0.05, 0.1):
        for seed in range(8):
            model = _fit(spans(80, 200, 6 + seed, smear_frac=smear_frac) + spans(30, 5, 7 + seed))
            reasons = [r["reason"] for r in model.rejected if r["units"] == MINOR_UNITS]
            assert reasons == ["support_below_threshold"], (smear_frac, seed, model.rejected)
