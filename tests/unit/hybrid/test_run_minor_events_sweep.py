"""Safety sweeps of the run-minority tier's within-run event test (heavy: full pipeline runs).

Both directions, on synthetic samples with realistic run stutter (``run_minor_synth``):

* a within-peak minority carrying insG, insG_pos58 or delinsAT (another base inside an X
  unit's C7 run) at or above the documented detection floor (``FLOOR_AF``,
  docs/reference/limitations.md) is never NEGATIVE, with two length peaks and with one;
  the HiFi-like two-peak insG minority at 30% of one allele with 2000 reads and at 40%
  with 600 reads (NEGATIVE before the within-run test) is included;
* wild-type samples whose reads carry G or AT noise inside random C runs, and
  homozygous normals with a site-specific G/AT-in-run artefact in 30-40% of reads (read
  for read a within-run event minority), are never PATHOGENIC.

The seeds and depths are the smallest that show each property; the 60-seed evidence
per claimed floor is in the calibration artefacts of the release validation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import run_minor_synth as rms
from tests.unit.hybrid import test_run_minor_events as tre
from tests.unit.hybrid import test_single_event as base

# Heavy synthetic safety sweep: its own CI job and make test-unit (never skipped).
pytestmark = pytest.mark.safety_sweep

S = DEFAULT_SETTINGS.hybrid
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
DEPTH = 2 * S.phase_max_site_reads
SEEDS = (0, 1)
# The shape found in review: HiFi-like two peaks, insG in one X unit, (depth, AF).
REVIEW_CELLS = ((S.phase_run_minor_max_reads, 0.30), (DEPTH, 0.40))
# Detection floor at DEPTH, the same for every within-run event, stutter shape and
# layout: 0 of 60 seeds NEGATIVE at it (docs/reference/limitations.md).
FLOOR_AF = 0.30
LAYOUTS = [(shape, layout) for shape in sorted(rms.SHAPES) for layout in rms.LAYOUTS]
ARTEFACT_SHARES = (0.30, 0.35, 0.40)
ARTEFACT_SHAPES = ("hifi", "ont")
ARTEFACT_DEPTH = S.phase_single_event_bound_reads
ARTEFACT_EVENTS = ("insG", "delinsAT")


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize(("depth", "af"), REVIEW_CELLS)
def test_review_shape_ins_g_minority_is_never_negative(
    tmp_path: Path, depth: int, af: float, seed: int
) -> None:
    records = rms.minority("hifi", af, depth, seed, "two_peak", "insG")
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize(("shape", "layout"), LAYOUTS)
@pytest.mark.parametrize("event", tre.IN_RUN_EVENTS)
def test_in_run_event_minority_at_the_floor_is_never_negative(
    tmp_path: Path, event: str, shape: str, layout: str, seed: int
) -> None:
    records = rms.minority(shape, FLOOR_AF, DEPTH, seed, layout, event)
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("layout", rms.LAYOUTS)
@pytest.mark.parametrize("shape", sorted(rms.SHAPES))
def test_wild_type_in_run_noise_is_never_pathogenic(
    tmp_path: Path, shape: str, layout: str, seed: int
) -> None:
    records = rms.in_run_noise(shape, tre.NOISE_RATE, DEPTH, seed, layout)
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != "PATHOGENIC", summary["hybrid"]


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("share", ARTEFACT_SHARES)
@pytest.mark.parametrize("shape", ARTEFACT_SHAPES)
@pytest.mark.parametrize("event", ARTEFACT_EVENTS)
def test_site_specific_in_run_artefact_is_never_pathogenic(
    tmp_path: Path, event: str, shape: str, share: float, seed: int
) -> None:
    records = rms.minority(shape, share, ARTEFACT_DEPTH, seed, "one_peak", event)
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != "PATHOGENIC", summary["hybrid"]["split_bases"]


# A seed of the artefact sweep above whose POA draft carried the artefact (C4 run): the
# single-event gate bounded only the non-draft (wild-type, 60%) side and the split
# called the 40% artefact group PATHOGENIC. Both sides are bounded now.
DRAFT_FOLLOWS_ARTEFACT = ("delinsAT", "hifi", 0.40, 7)


def test_draft_following_an_artefact_is_never_pathogenic(tmp_path: Path) -> None:
    event, shape, share, seed = DRAFT_FOLLOWS_ARTEFACT
    records = rms.minority(shape, share, ARTEFACT_DEPTH, seed, "one_peak", event)
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != "PATHOGENIC", summary["hybrid"]["split_bases"]


# Seeds of the saturating ONT-like artefact sweep: wild-type reads misaligned next to
# the stuttering run showed neither allele at the delinsAT column and were skipped, so
# both shares of a 40% artefact were bounded above the floor. Every observing read
# counts now.
NEITHER_ALLELE_SEEDS = (20, 29)


@pytest.mark.parametrize("seed", NEITHER_ALLELE_SEEDS)
def test_reads_with_neither_allele_do_not_inflate_an_artefact(tmp_path: Path, seed: int) -> None:
    records = rms.minority("ont_saturating", 0.40, ARTEFACT_DEPTH, seed, "one_peak", "delinsAT")
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != "PATHOGENIC", summary["hybrid"]["split_bases"]
