"""Task 15h: PCR dimer peaks and inter-allele smear are length artefacts, not alleles.

Safety first: a real allele is never called a dimer (no read of a real allele carries a
motif-9 -> motif-1 amplicon junction), accepted peaks are never touched, a dimer
explanation that does not fit (too many dimer reads for a minority artefact, or a
residual cluster left after removing the dimer reads) keeps the sample unresolved, and a
real minor peak between two alleles is not silently called smear.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from muc_one_span.hybrid import lengths as lengths_module
from muc_one_span.hybrid.dimers import DimerCall, find_junction, recognise_dimers
from muc_one_span.hybrid.lengths import GATE_RELEVANT_REJECTIONS, LengthModel, fit_length_model
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S)
UNIT = ANCH.unit_bp
ERR = 0.02
SHORT_INNER = 30  # 39 units in total (synth.PRE + inner + synth.POST)
LONG_INNER = 45
SHORT_UNITS = len(synth.PRE) + SHORT_INNER + len(synth.POST)
LONG_UNITS = len(synth.PRE) + LONG_INNER + len(synth.POST)
MAJOR_READS = 400
DIMER_READS = 8  # 2 % of MAJOR_READS: a minority artefact at the default parent fraction


def _allele(inner: int) -> str:
    return synth.allele(["X"] * inner)


def _spans(inner: int, n: int, seed: int, **kw: float) -> list[SpanRead]:
    reads = synth.reads(_allele(inner), n, err=ERR, seed=seed, **kw)
    return categorize_reads(reads, ANCH, S).spanning


def _dimer_spans(inner_a: int, inner_b: int, n: int, seed: int) -> list[SpanRead]:
    """``n`` spanning dimer products (motif 1 of copy a .. motif 9 of copy b)."""
    floor = len(_allele(inner_a)) + len(_allele(inner_b))
    out: list[SpanRead] = []
    batch = 0
    while len(out) < n:
        reads = synth.concatemers(
            _allele(inner_a), _allele(inner_b), 4 * n, err=ERR, seed=seed + batch
        )
        out += [s for s in categorize_reads(reads, ANCH, S).spanning if s.length >= floor]
        batch += 1
    return out[:n]


def _fit(spans: list[SpanRead], settings: HybridSettings = S) -> LengthModel:
    return fit_length_model(spans, settings, ANCH)


def _units(model: LengthModel) -> list[int]:
    return sorted(round(p.center_bp / UNIT) for p in model.peaks)


def _dimers(model: LengthModel) -> list[dict[str, object]]:
    return [r for r in model.rejected if r["reason"] == "dimer"]


# --- junction detection -------------------------------------------------------------


def test_junction_found_in_dimer_reads_with_both_parts() -> None:
    dimers = _dimer_spans(SHORT_INNER, LONG_INNER, 10, seed=1)
    for span in dimers:
        junction = find_junction(span, ANCH, S)
        assert junction is not None
        assert round(junction.first_bp / UNIT) == SHORT_UNITS
        assert round(junction.second_bp / UNIT) == LONG_UNITS


def test_junction_never_found_in_genuine_allele_reads() -> None:
    for inner, seed in ((SHORT_INNER, 2), (LONG_INNER, 3), (2 * SHORT_INNER + 9, 4)):
        genuine = _spans(inner, 60, seed, smear_frac=0.3)
        assert [find_junction(s, ANCH, S) for s in genuine] == [None] * len(genuine)


# --- dimer recognition --------------------------------------------------------------


def test_homozygous_dimer_peak_is_an_artefact_not_a_rejected_peak() -> None:
    major = _spans(SHORT_INNER, MAJOR_READS, 5)
    dimers = _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 6)
    before = _fit(major + dimers, replace(S, dimer_recognition=False))
    assert before.gate_relevant_rejections != []  # today's INCONCLUSIVE driver
    model = _fit(major + dimers)
    assert _units(model) == [SHORT_UNITS]
    assert model.gate_relevant_rejections == []
    [entry] = _dimers(model)
    assert entry["support"] == DIMER_READS
    assert entry["parent_units"] == [SHORT_UNITS, SHORT_UNITS]
    assert sorted(s.name for s in model.dimer_products) == sorted(s.name for s in dimers)
    members = {m.name for p in model.peaks for m in p.members}
    assert members.isdisjoint(s.name for s in dimers)  # excluded from consensus/evidence
    assert not {s.name for s in model.unassigned} & {s.name for s in dimers}


def test_heterozygous_twice_and_sum_dimers_are_recognised() -> None:
    reads = (
        _spans(SHORT_INNER, MAJOR_READS, 7)
        + _spans(LONG_INNER, MAJOR_READS, 8)
        + _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 9)
        + _dimer_spans(SHORT_INNER, LONG_INNER, DIMER_READS, 10)
        + _dimer_spans(LONG_INNER, LONG_INNER, DIMER_READS, 11)
    )
    model = _fit(reads)
    assert _units(model) == [SHORT_UNITS, LONG_UNITS]
    assert model.gate_relevant_rejections == []
    parents = sorted(tuple(r["parent_units"]) for r in _dimers(model))  # type: ignore[arg-type]
    assert parents == [
        (SHORT_UNITS, SHORT_UNITS),
        (SHORT_UNITS, LONG_UNITS),
        (LONG_UNITS, LONG_UNITS),
    ]


def test_weak_real_allele_at_twice_the_length_is_never_a_dimer() -> None:
    # A genuine long allele exactly (and near) 2 x L with PCR dropout: its reads carry no
    # amplicon junction, so it stays a gate-relevant rejected peak (INCONCLUSIVE).
    for extra_units in (0, 1, -1):
        inner = SHORT_UNITS * 2 - len(synth.PRE) - len(synth.POST) + extra_units
        model = _fit(_spans(SHORT_INNER, MAJOR_READS, 12) + _spans(inner, DIMER_READS, 13))
        assert _dimers(model) == []
        assert [r["reason"] for r in model.gate_relevant_rejections] == ["support_below_threshold"]


def test_weak_real_allele_mixed_with_dimers_at_twice_the_length_stays_unresolved() -> None:
    inner = SHORT_UNITS * 2 - len(synth.PRE) - len(synth.POST)
    reads = (
        _spans(SHORT_INNER, MAJOR_READS, 14)
        + _spans(inner, DIMER_READS, 15)
        + _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 16)
    )
    model = _fit(reads)
    assert _units(model) == [SHORT_UNITS]
    assert [(r["reason"], r["support"]) for r in model.gate_relevant_rejections] == [
        ("support_below_threshold", DIMER_READS)
    ]


def test_accepted_allele_at_twice_the_length_is_kept_with_its_dimers_removed() -> None:
    inner = SHORT_UNITS * 2 - len(synth.PRE) - len(synth.POST)
    reads = (
        _spans(SHORT_INNER, MAJOR_READS, 17)
        + _spans(inner, MAJOR_READS // 2, 18)
        + _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 19)
    )
    model = _fit(reads)
    assert _units(model) == [SHORT_UNITS, 2 * SHORT_UNITS]
    assert model.gate_relevant_rejections == []


def test_peak_made_of_dimer_reads_is_not_accepted_as_an_allele() -> None:
    # Enough dimer reads to clear the far-peak threshold (today: a false second allele at
    # 2 x L), still a minority artefact relative to the parent.
    n_dimers = 2 * DIMER_READS
    major = _spans(SHORT_INNER, MAJOR_READS, 28)
    dimers = _dimer_spans(SHORT_INNER, SHORT_INNER, n_dimers, 29)
    assert len(_fit(major + dimers, replace(S, dimer_recognition=False)).peaks) == 2
    model = _fit(major + dimers)
    assert _units(model) == [SHORT_UNITS]
    assert model.gate_relevant_rejections == []
    assert [r["support"] for r in _dimers(model)] == [n_dimers]


def test_dimer_reads_above_the_parent_fraction_are_not_explained() -> None:
    # More dimer reads than a minority artefact allows: ambiguous -> unresolved.
    major = _spans(SHORT_INNER, MAJOR_READS, 20)
    dimers = _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 21)
    parent = next(p.support for p in _fit(major).peaks)
    strict = replace(S, dimer_max_parent_frac=(DIMER_READS - 1) / parent)
    model = _fit(major + dimers, strict)
    assert _dimers(model) == []
    assert [r["reason"] for r in model.gate_relevant_rejections] == ["support_below_threshold"]
    assert model.dimer_products == []


def test_dimer_recognition_can_be_switched_off() -> None:
    major = _spans(SHORT_INNER, MAJOR_READS, 22)
    dimers = _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 23)
    model = _fit(major + dimers, replace(S, dimer_recognition=False))
    assert _dimers(model) == [] and model.dimer_products == []
    assert model.gate_relevant_rejections != []


def test_three_real_alleles_with_dimers_stay_max_alleles() -> None:
    reads = (
        _spans(SHORT_INNER, 200, 24)
        + _spans(LONG_INNER, 200, 25)
        + _spans(LONG_INNER + 15, 200, 26)
        + _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 27)
    )
    model = _fit(reads)
    assert "max_alleles" in [r["reason"] for r in model.gate_relevant_rejections]


# --- inter-allele smear -------------------------------------------------------------

INTER_LONG_INNER = 80
SMEAR_SEEDS = 10


def _inter(seed: int, **extra: list[SpanRead]) -> LengthModel:
    reads = _spans(SHORT_INNER, MAJOR_READS, seed) + _spans(
        INTER_LONG_INNER, 250, seed + 50, smear_frac=0.3
    )
    return _fit(reads + extra.get("minor", []))


def test_smear_between_the_alleles_is_smear_not_a_rejected_peak() -> None:
    for seed in range(SMEAR_SEEDS):
        model = _inter(seed)
        assert len(model.peaks) == 2, (seed, model.rejected)
        assert model.gate_relevant_rejections == [], (seed, model.rejected)
        inter = [r for r in model.rejected if r.get("smear_region") == "inter_allele"]
        assert inter and {r["reason"] for r in inter} == {"smear"}


def test_inter_allele_smear_test_can_be_switched_off() -> None:
    reads = _spans(SHORT_INNER, MAJOR_READS, 0) + _spans(INTER_LONG_INNER, 250, 50, smear_frac=0.3)
    model = _fit(reads, replace(S, smear_test_inter_allele=False))
    assert model.gate_relevant_rejections != []


def test_real_minor_between_the_alleles_is_never_silent_smear() -> None:
    # A contamination/mosaic-like third peak between the alleles (out of scope to call)
    # must stay gate-relevant: significant excess over the smear -> INCONCLUSIVE.
    minor_units = len(synth.PRE) + 55 + len(synth.POST)
    for minor_n in (15, 40, 120):
        for seed in range(SMEAR_SEEDS):
            minor = _spans(55, minor_n, seed + 900)
            model = _inter(seed, minor=minor)
            flagged = [
                r
                for r in model.rejected
                if r["units"] == minor_units and r["reason"] in GATE_RELEVANT_REJECTIONS
            ]
            assert flagged, (minor_n, seed, model.rejected)


# --- fail-closed edges --------------------------------------------------------------


def test_no_accepted_peak_means_no_dimer() -> None:
    dimers = _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 30)
    assert recognise_dimers(dimers, [], ANCH, S) == DimerCall([], [], [])


def test_dimer_of_an_allele_that_was_not_accepted_is_not_explained() -> None:
    # Junction parts must both match accepted peaks: a (short, unseen) dimer stays.
    major = _spans(SHORT_INNER, MAJOR_READS, 31)
    foreign = _dimer_spans(SHORT_INNER, LONG_INNER + 30, DIMER_READS, 32)
    model = _fit(major + foreign)
    assert _dimers(model) == [] and model.dimer_products == []
    assert model.gate_relevant_rejections != []


def test_refit_that_loses_a_parent_keeps_the_first_fit(monkeypatch: pytest.MonkeyPatch) -> None:
    major = _spans(SHORT_INNER, MAJOR_READS, 33)
    dimers = _dimer_spans(SHORT_INNER, SHORT_INNER, DIMER_READS, 34)
    off = _fit(major + dimers, replace(S, dimer_recognition=False))
    far_away = float(10 * MAJOR_READS * UNIT)  # a parent no refit peak can be near

    def fake(*_: object) -> DimerCall:
        return DimerCall(dimers, [{"reason": "dimer"}], [far_away])

    monkeypatch.setattr(lengths_module, "recognise_dimers", fake)
    model = _fit(major + dimers)
    assert model.rejected == off.rejected and model.dimer_products == []
