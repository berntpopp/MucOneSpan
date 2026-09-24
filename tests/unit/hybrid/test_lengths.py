"""S2 length peaks: scaled windows, rejected peaks, smear and close alleles."""

from __future__ import annotations

from muc_one_span.hybrid.lengths import fit_length_model
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD)


def spans(inner_units: int, n: int, seed: int, **kw: float) -> list[SpanRead]:
    seq = synth.allele(["X"] * inner_units)
    return categorize_reads(synth.reads(seq, n, err=0.02, seed=seed, **kw), ANCH, S).spanning


def test_two_distant_alleles_with_minor_long_allele() -> None:
    model = fit_length_model(spans(30, 200, 1) + spans(70, 12, 2), S)
    assert sorted(round(p.center_bp / 60) for p in model.peaks) == [39, 79]
    assert model.gate_relevant_rejections == []


def test_smear_is_short_product_not_allele() -> None:
    model = fit_length_model(spans(60, 120, 3, smear_frac=0.45), S)
    assert len(model.peaks) == 1
    assert model.short_product_fraction > 0.3
    assert model.gate_relevant_rejections == []


def test_real_short_allele_amid_smear_is_kept() -> None:
    model = fit_length_model(spans(30, 40, 8) + spans(60, 120, 9, smear_frac=0.3), S)
    assert sorted(round(p.center_bp / 60) for p in model.peaks) == [39, 69]


def test_close_alleles_are_not_silently_dropped() -> None:
    # Alleles one unit apart must yield two peaks or a large unassigned fraction.
    model = fit_length_model(spans(40, 150, 4) + spans(41, 150, 5), S)
    assert len(model.peaks) == 2 or model.unassigned_fraction > 0.2


def test_rejected_minor_peak_is_reported() -> None:
    model = fit_length_model(spans(30, 200, 6) + spans(80, 5, 7), S)
    assert len(model.peaks) == 1
    assert [(r["support"], r["reason"]) for r in model.gate_relevant_rejections] == [
        (5, "support_below_threshold")
    ]


def test_third_real_peak_is_max_alleles() -> None:
    model = fit_length_model(spans(30, 60, 10) + spans(50, 60, 11) + spans(70, 60, 12), S)
    assert len(model.peaks) == 2
    assert [r["reason"] for r in model.gate_relevant_rejections] == ["max_alleles"]


def test_single_outlier_read_is_noise_not_gate_relevant() -> None:
    model = fit_length_model(spans(30, 100, 13) + spans(90, 1, 14), S)
    assert [r["reason"] for r in model.rejected] == ["noise"]
    assert model.gate_relevant_rejections == []
