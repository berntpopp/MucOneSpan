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


def test_smear_never_silently_becomes_a_second_peak() -> None:
    # I1 (fix round 1) property test: a single real allele plus PCR/sequencing smear must
    # never yield a second *accepted* peak with no gate-relevant signal anywhere in the
    # model. Sweeps the brief's own scenario (spans(60, 120, seed, smear_frac=...)) across
    # seeds 0-199 at both smear fractions named in the controller's ruling.
    failures = []
    for smear_frac in (0.3, 0.45):
        for seed in range(200):
            model = _fit(spans(60, 120, seed, smear_frac=smear_frac))
            if len(model.peaks) == 2 and model.gate_relevant_rejections == []:
                failures.append((smear_frac, seed))
    assert failures == []
