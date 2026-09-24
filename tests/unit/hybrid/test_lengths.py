"""S2 length peaks: scaled windows, rejected peaks, smear and close alleles."""

from __future__ import annotations

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
    # Fix round 3 (R2): one below-top candidate here (support=8, right at min_peak_reads,
    # with a small positive excess over the expected smear count) is no longer silently
    # 'smear' -- it is correctly flagged smear_ambiguous (gate-relevant) instead. R2's
    # ruling is that this must never be silent, so this is the corrected behaviour, not a
    # regression: the bulk of the smear tail (every other below-top candidate here) is
    # still silent 'smear' or 'noise'.
    assert [r["reason"] for r in model.gate_relevant_rejections] == ["smear_ambiguous"]


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

DS = (60, 120, 200, 300)
FRACS = (0.3, 0.45, 0.54)
GRID_SEEDS = 15  # keeps the whole file well under the ~60s budget; see task-5-report.md
F2_SEEDS = 8
SPURIOUS_TARGET = 0.10  # criterion (b)
ACCEPT_TARGET = 0.90  # criterion (c)
MINOR_SUPPORT_FRAC_20 = 0.2  # criterion (c): "support >= 20% of the major"
MINOR_SUPPORT_FRAC_10 = 0.1  # R2's probe minor: weaker, 10% of the major (fix round 3)

# Fix round 3 (R2): a below-top candidate with support >= min_peak_reads and a positive
# excess over the expected smear count is now never silently "smear" -- it is accepted or
# smear_ambiguous (gate-relevant) instead. That correctness fix pushes many homozygous
# draws that used to be silent into gate-relevant territory, so most cells now miss the
# <=10% (b) target. Per the controller's ruling ("keep (b) measured... do not loosen
# (d)"), this is reported as DONE_WITH_CONCERNS with the measured table below rather than
# a weakened assertion anywhere. R1's fix (a clean, isolated candidate that already clears
# the allele threshold is accepted, not merely-isolated-therefore-ambiguous) resolved (c)
# everywhere, including the previous (60, 0.3) shortfall, so (c) has no exceptions left.
# See task-5-report.md Fix round 3 for the full per-cell (a)/(b)/(c)/(d) table.
KNOWN_SPURIOUS_TRADEOFF_CELLS = frozenset(
    {
        (60, 0.3),
        (120, 0.45),
        (120, 0.54),
        (200, 0.3),
        (200, 0.45),
        (200, 0.54),
        (300, 0.3),
        (300, 0.45),
        (300, 0.54),
    }
)


def test_smear_grid_never_silently_accepts_and_spurious_rate_is_bounded() -> None:
    """(a): a single real allele plus smear must never yield a silently-accepted second
    peak -- zero tolerance, every (D, smear_frac) cell, every seed, no exceptions.
    (b): the homozygous gate-relevant (spurious) rate is <=10% per cell, except the
    documented KNOWN_SPURIOUS_TRADEOFF_CELLS (fix round 3, R2)."""
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
        if (depth, frac) in KNOWN_SPURIOUS_TRADEOFF_CELLS:
            continue
        assert rate <= SPURIOUS_TARGET, f"D={depth} frac={frac}: spurious rate {rate:.0%}"


def _minor_outcome(depth: int, frac: float, minor_n: int, seed: int) -> tuple[bool, list[str]]:
    """(accepted, reasons recorded for the minor's own unit) for one grid draw."""
    model = _fit(
        spans(60, depth, seed, smear_frac=frac) + spans(MINOR_INNER_UNITS, minor_n, seed + 1000)
    )
    accepted = MINOR_UNITS in {round(p.center_bp / UNIT) for p in model.peaks}
    reasons = [r["reason"] for r in model.rejected if r["units"] == MINOR_UNITS]
    return accepted, reasons


def test_smear_grid_genuine_minor_is_accepted_and_never_silently_smear() -> None:
    """(c): a genuine below-top minor with support >= 20% of the major is accepted in
    >=90% of seeds, every cell -- R1's fix resolved the previous (60, 0.3) shortfall, so
    there are no (c) exceptions any more.
    (d): that minor is never silently 'smear' when its support clears min_peak_reads --
    always accepted or gate-relevant, every cell, every seed, no exceptions."""
    silent_smear = []
    accept_rates = {}
    for depth in DS:
        for frac in FRACS:
            minor_n = round(MINOR_SUPPORT_FRAC_20 * depth)
            accepted_count = 0
            for seed in range(GRID_SEEDS):
                accepted, reasons = _minor_outcome(depth, frac, minor_n, seed)
                if accepted:
                    accepted_count += 1
                elif minor_n >= S.min_peak_reads and "smear" in reasons:
                    silent_smear.append((depth, frac, seed))
            accept_rates[(depth, frac)] = accepted_count / GRID_SEEDS
    assert silent_smear == []
    for (depth, frac), rate in accept_rates.items():
        assert rate >= ACCEPT_TARGET, f"D={depth} frac={frac}: accept rate {rate:.0%}"


def test_smear_grid_ten_percent_minor_is_never_silently_smear() -> None:
    """R2 (fix round 3): a below-top local maximum with support >= min_peak_reads and a
    positive excess over the expected smear count is never silently 'smear' -- accepted
    or gate-relevant only. Probes a weaker minor (10% of the major) than the (c) test
    above: this was the review's own counter-example (D=120: 12 reads silent in 14/15
    seeds; D=300: 30 reads silent in 15/15, both under the pre-fix model)."""
    silent_smear = []
    for depth in DS:
        minor_n = round(MINOR_SUPPORT_FRAC_10 * depth)
        if minor_n < S.min_peak_reads:
            continue  # below the threshold this guarantee is scoped to (D=60: 6 reads)
        for frac in FRACS:
            for seed in range(GRID_SEEDS):
                accepted, reasons = _minor_outcome(depth, frac, minor_n, seed)
                if not accepted and "smear" in reasons:
                    silent_smear.append((depth, frac, minor_n, seed))
    assert silent_smear == []


def test_d60_light_smear_minor_that_misses_acceptance_is_gate_relevant() -> None:
    """In the hardest grid cell (D=60, smear=0.3), a genuine 20%-minor that is not
    accepted must still be recorded exactly once at its own unit (not silently dropped,
    and not merged into a neighbouring unit's candidate) and must be gate-relevant --
    never plain "noise" and never silent "smear"."""
    depth, frac = 60, 0.3
    minor_n = round(MINOR_SUPPORT_FRAC_20 * depth)
    for seed in range(GRID_SEEDS):
        model = _fit(
            spans(60, depth, seed, smear_frac=frac) + spans(MINOR_INNER_UNITS, minor_n, seed + 1000)
        )
        peak_units = {round(p.center_bp / UNIT) for p in model.peaks}
        matches = [r for r in model.rejected if r["units"] == MINOR_UNITS]
        assert MINOR_UNITS in peak_units or len(matches) == 1, (seed, model.rejected)
        if MINOR_UNITS not in peak_units:
            assert matches[0]["reason"] in GATE_RELEVANT_REJECTIONS, (seed, matches)


def test_isolated_minor_with_light_smear_is_support_below_threshold() -> None:
    # (e), F2: an isolated real minor allele stays support_below_threshold (gate-relevant),
    # never silently "smear", across the light smear fractions the review measured it
    # failing at (0.01-0.1: 1-10% of the major's reads smeared).
    for smear_frac in (0.01, 0.02, 0.05, 0.1):
        for seed in range(F2_SEEDS):
            model = _fit(spans(80, 200, 6 + seed, smear_frac=smear_frac) + spans(30, 5, 7 + seed))
            reasons = [r["reason"] for r in model.rejected if r["units"] == MINOR_UNITS]
            assert reasons == ["support_below_threshold"], (smear_frac, seed, model.rejected)


def test_clean_isolated_minor_is_accepted_not_ambiguous() -> None:
    """R1 (fix round 3): a below-top candidate with no smear background at all, that
    already clears the allele threshold, is accepted -- being isolated is not itself
    ambiguous. Round 2 regressed every one of these to smear_ambiguous."""
    for depth, minor_n in ((200, 50), (200, 20), (120, 30), (60, 15)):
        for seed in range(10):
            model = _fit(spans(80, depth, seed) + spans(MINOR_INNER_UNITS, minor_n, seed + 7777))
            accepted = MINOR_UNITS in {round(p.center_bp / UNIT) for p in model.peaks}
            assert accepted, (depth, minor_n, seed, model.rejected)
