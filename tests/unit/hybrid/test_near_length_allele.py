"""Task 15i: the smear relabelling must never absorb a peak that is a real allele.

Reproduces the shape of a simulated ONT amplicon case (a regression shape from the
unsealed test split, used non-blind): alleles of 96 and 93 repeat units with almost the
same support (537 and 531 reads) over heavy PCR smear. The 93-unit peak was the
shorter of the two and ~2.9 units below the top, so it was smear-tested in the
below-top region. That region ends ``smear_short_product_units`` below the top, just
past the candidate's own assignment window, so the candidate's right background side
was clipped to a sliver (0.36 bp, no reads). With 441 core reads, "same density as a
0.36 bp side holding no reads" still has p = 0.025 > alpha: the sliver, which cannot
tell smear from an allele, decided the verdict, the real allele became silent 'smear'
and the sample was NEGATIVE with the variant-carrying allele lost.

Fixes: (1) a background side clipped by the region edge that holds fewer than
``smear_background_min_reads`` reads and could not make the candidate significant even
with no reads carries no information and is skipped while the other side can decide; (2) a candidate with at least
``smear_guard_top_frac`` of the top peak's reads is never smear-tested (it faces the
allele support rules: accepted, or gate-relevant). A guard at the allele threshold
(``max(min_peak_reads, frac * N)``) was evaluated and not adopted: at low depth that
threshold is a handful of reads and smear debris clusters reach it, so homozygous
smear would become gate-relevant far beyond the Task 5 criterion (b).
"""

from __future__ import annotations

import dataclasses
import random

import pytest

from muc_one_span.hybrid import lengths
from muc_one_span.hybrid.lengths import fit_length_model
from muc_one_span.hybrid.smear import smear_test, smear_verdict
from muc_one_span.hybrid.spans import Anchors, SpanRead
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid.synth import LAYOUT

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S, LAYOUT)
UNIT = ANCH.unit_bp
# The case: allele peaks (bp, reads) and its smear (reads below the shorter allele).
# The case's KDE centres sat 20 and 16 bp short of 96 and 93 x unit (176 bp apart).
TOP = (96 * UNIT - 20, 537)
NEAR = (93 * UNIT - 16, 531)
SMEAR_READS = 790
SMEAR_GAP_BP = 264  # the case's smear ended ~4.4 units below the 93-unit allele
SPAN_SD_BP = 15.0  # ONT span-length noise of the case's peaks (about +-40 bp at 2 sd)
SEEDS = (0, 1, 2)
# Shifts of the shorter allele around the case geometry; at 0 the sliver side decided
# the verdict for seeds 0 and 2 before the fix.
OFFSETS_BP = (-8, -4, 0, 4, 8)


def _span(i: int, length: int) -> SpanRead:
    return SpanRead(f"s{i}", "A" * length, 20.0, "+", 0, "motif")


def _case_lengths(seed: int, near_reads: int = NEAR[1], offset: int = 0) -> list[float]:
    rng = random.Random(seed)
    out = [rng.gauss(TOP[0], SPAN_SD_BP) for _ in range(TOP[1])]
    out += [rng.gauss(NEAR[0] + offset, SPAN_SD_BP) for _ in range(near_reads)]
    lo = S.min_span_units * UNIT
    out += [rng.uniform(lo, NEAR[0] - SMEAR_GAP_BP) for _ in range(SMEAR_READS)]
    return [float(round(x)) for x in out]


def _model(lens: list[float], settings: HybridSettings = S) -> lengths.LengthModel:
    return fit_length_model([_span(i, int(x)) for i, x in enumerate(lens)], settings, ANCH)


@pytest.mark.parametrize("offset", OFFSETS_BP)
@pytest.mark.parametrize("seed", SEEDS)
def test_near_length_allele_below_the_top_is_kept(seed: int, offset: int) -> None:
    """The regression: a real 93-unit allele 3 units below a 96-unit top is an allele."""
    model = _model(_case_lengths(seed, offset=offset))
    units = sorted(round(p.center_bp / UNIT) for p in model.peaks)
    assert len(units) == 2 and abs(units[0] - 93) <= 1 and units[1] == 96, model.rejected[:3]


def test_case_geometry_reproduces_the_uninformative_sliver_side() -> None:
    """Why the smear test alone fails here (documents the residual, see limitations):
    the right side is clipped to a sliver with no reads, which cannot reach significance
    but still sets the p value above alpha."""
    lens = _case_lengths(SEEDS[0])
    sel = lengths._select(lens, S, UNIT)
    top = max(sel.support, key=lambda c: sel.support[c])
    near = next(c for c in sel.support if abs(c / UNIT - 93) < 1)
    region = (min(min(lens), S.min_span_units * UNIT), top - S.smear_short_product_units * UNIT)
    test = smear_test(near, lengths.window_bp(near, S, UNIT), lens, region, S, UNIT)
    sliver = min(test.sides, key=lambda side: side.width_bp)
    assert sliver.reads == 0 and sliver.width_bp < SPAN_SD_BP
    assert smear_verdict(test.p_value, 1, S) != "significant"


def test_without_the_guard_the_case_allele_is_lost() -> None:
    """The guard is what keeps it: at smear_guard_top_frac=1 the case is smear again."""
    unguarded = dataclasses.replace(S, smear_guard_top_frac=1.0)
    model = _model(_case_lengths(SEEDS[0]), unguarded)
    assert len(model.peaks) == 1
    assert any(abs(r["units"] - 93) <= 1 and r["reason"] == "smear" for r in model.rejected)


def _always(verdict: str, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(lengths, "smear_verdict", lambda *_a, **_k: verdict)


def test_candidate_with_top_level_support_is_never_smear(monkeypatch: pytest.MonkeyPatch) -> None:
    """Whatever the smear test says, >= smear_guard_top_frac of the top is an allele."""
    _always("smear", monkeypatch)
    model = _model(_case_lengths(SEEDS[0]))
    assert len(model.peaks) == 2
    assert all(
        r["reason"] != "smear" or r["support"] < S.smear_guard_top_frac * TOP[1]
        for r in model.rejected
    )


def test_candidate_below_the_guard_is_left_to_the_smear_test(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Below the guard the smear test decides, whatever the absolute read count (Task 5
    ruling: smear debris grows with depth, so a read count alone never makes an allele)."""
    _always("smear", monkeypatch)
    near_reads = int(S.smear_guard_top_frac * TOP[1]) - 1
    model = _model(_case_lengths(SEEDS[0], near_reads))
    near = [r for r in model.rejected if abs(r["units"] - 93) <= 1]
    assert near and near[0]["reason"] == "smear", model.rejected[:3]


def test_inter_allele_relabel_respects_the_guard(monkeypatch: pytest.MonkeyPatch) -> None:
    """The inter-allele smear test never relabels a candidate with top-level support."""
    sel = lengths._Selection(
        kept=[1800.0, 4200.0],
        support={1800.0: 400, 4200.0: 300, 3000.0: 200},
        rejected=[(3000.0, {"units": 50, "support": 200, "reason": "support_below_threshold"})],
    )
    _always("smear", monkeypatch)
    lengths._inter_allele_smear(sel, [1800.0] * 400 + [4200.0] * 300 + [3000.0] * 200, S, UNIT)
    assert sel.rejected[0][1]["reason"] == "support_below_threshold"


def test_guard_setting_is_validated() -> None:
    for bad in (0, -0.1, 1.5):
        with pytest.raises(ValueError, match="smear_guard_top_frac"):
            HybridSettings(smear_guard_top_frac=bad)
