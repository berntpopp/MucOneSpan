"""Task 15l safety sweeps (heavy: full pipeline runs; grouped here for a dedicated CI job).

Both directions of the run-minority tier (``run_minor``) and the single-event share
floor (``phase_single_event_min_share``), on synthetic samples with realistic run
stutter (``run_minor_synth``: HiFi-like length-dependent, ONT-like strand-asymmetric and
ONT-like saturating stutter; poor reads stutter twice as often):

* a within-peak dupC minority at or above the documented tier-only detection floor
  (``FLOOR_AF``, docs/reference/limitations.md) is never NEGATIVE, with two length
  peaks and with one (the 60-seed evidence is in the Task 15l report);
* wild-type samples under the same stutter are never PATHOGENIC, and the tier never
  turns one of them from NEGATIVE into INCONCLUSIVE;
* a homozygous normal whose C unit's C6 run carries a site-specific +1 C artefact in
  30-40% of reads (the simulated HiFi ``H1_hifi`` shape; the artefact reads are exactly
  dupA carrier reads) is never PATHOGENIC;
* the 15k two-peak all-low-quality dupC minority shapes (15% of one allele) are never
  NEGATIVE (``test_quality_known_events``; no longer expected failures).

The depths, seeds and read counts are the smallest that show each property; the full
grids (AF 0.10-0.30, depth 60-2000, five seeds) are in the Task 15l report.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import run_minor_synth as rms
from tests.unit.hybrid import test_single_event as base

S = DEFAULT_SETTINGS.hybrid
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
# The event peak fills the phase sample (phase_max_site_reads reads per allele).
DEPTH = 2 * S.phase_max_site_reads
SEEDS = (0, 1)
# Tier-only detection floor at DEPTH (monotone envelope over depth >= DEPTH), per
# (stutter shape, layout): the run-minority tier itself fired and nothing was NEGATIVE
# in 60 seeds (Task 15l fix round 1; docs/reference/limitations.md). ONT-like
# strand-asymmetric single peaks have no tier floor up to 0.40 (their minorities are
# blocked by other gates only) and are not hard-tested here.
FLOOR_AF = {
    ("hifi", "two_peak"): 0.50,
    ("hifi", "one_peak"): 0.35,
    ("ont", "two_peak"): 0.30,
    ("ont_saturating", "two_peak"): 0.25,
    ("ont_saturating", "one_peak"): 0.20,
}
# Share of reads carrying the wild-type run artefact (the H1_hifi shape: 0.327).
ARTEFACT_SHARES = (0.30, 0.35, 0.40)
ARTEFACT_SHAPES = ("hifi", "ont")
# The artefact depth gives the single-event share bound its full sample.
ARTEFACT_DEPTH = S.phase_single_event_bound_reads
TIER_OFF = dataclasses.replace(
    DEFAULT_SETTINGS, hybrid=dataclasses.replace(S, phase_run_minor_scope="off")
)


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize(("shape", "layout"), sorted(FLOOR_AF))
def test_run_minority_at_the_floor_is_never_negative(
    tmp_path: Path, shape: str, layout: str, seed: int
) -> None:
    records = rms.minority(shape, FLOOR_AF[(shape, layout)], DEPTH, seed, layout)
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("layout", rms.LAYOUTS)
@pytest.mark.parametrize("shape", sorted(rms.SHAPES))
def test_wild_type_under_realistic_stutter_is_never_pathogenic(
    tmp_path: Path, shape: str, layout: str, seed: int
) -> None:
    records = rms.minority(shape, 0.0, DEPTH, seed, layout)
    for sub in ("on", "off"):
        (tmp_path / sub).mkdir()
    summary, decision = base._run(tmp_path / "on", records)
    assert decision["state"] != "PATHOGENIC", summary["hybrid"]
    if "unconfirmed_run_minor" in summary["hybrid"]["split_bases"]:
        # The tier costs no wild-type NEGATIVE: such a sample is INCONCLUSIVE without it.
        _summary, off = base._run(tmp_path / "off", records, TIER_OFF)
        assert off["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("share", ARTEFACT_SHARES)
@pytest.mark.parametrize("shape", ARTEFACT_SHAPES)
def test_wild_type_run_artefact_is_never_pathogenic(
    tmp_path: Path, shape: str, share: float, seed: int
) -> None:
    summary, decision = base._run(tmp_path, rms.artefact(shape, share, ARTEFACT_DEPTH, seed))
    assert decision["state"] != "PATHOGENIC", summary["hybrid"]["split_bases"]
