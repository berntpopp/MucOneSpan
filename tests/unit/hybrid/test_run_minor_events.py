"""A within-run event minority (insG, insG_pos58, delinsAT) blocks NEGATIVE like dupC.

insG, insG_pos58 and delinsAT put another base inside an X unit's C7 run
(``CCCCCCGC``, ``CCCCCGCC``, ``CCATCCCC``). A carrier read never observes that run
cleanly (``polish.run_observation``), so the run-minority tier, which reads clean
observations only, could not see such a minority: a HiFi-like two-peak insG minority at
30% of one allele with 2000 reads, or 40% with 600, was NEGATIVE. The site table sees
the carrier read as the run shortened to its longest C stretch (C7 -> C6, C5 or C4),
the template's run signature (``known_events``).

The tier therefore also tests each run's bounded but impure observations (both
bounding bases kept, another base inside the run): the share of reads whose longest
stretch is a known-event run signature, against the rate at the peer runs, is a
two-component mixture whose minor component reads the signature length with the
run-length stutter probability of that length. Its one-sided lower bound must reach the
same ``phase_run_minor_min_share``.
"""

from __future__ import annotations

import dataclasses
import random
from typing import Any

import pytest

from muc_one_span.hybrid.known_events import known_event_sites, template_signatures
from muc_one_span.hybrid.phase_sites import features
from muc_one_span.hybrid.run_minor import run_minor_sites
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import run_minor_synth as rms
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_run_minor as trm
from tests.unit.hybrid import test_single_event as base

S = DEFAULT_SETTINGS.hybrid
KNOWN = known_event_sites(synth.RD, S)
IN_RUN_EVENTS = ("insG", "insG_pos58", "delinsAT")
SEED = 0
# In-run noise in this share of reads (one random C run each) is ordinary sequencing
# noise: it is spread over every run, so the peer runs show it too.
NOISE_RATE = 0.3


def _signature(event: str) -> tuple[str, int, int]:
    """The (base, parent length, carrier length) run signature of ``event`` in X."""
    context = synth.RD.repeats[synth.RD.canonical_repeat]
    runs, _ins = template_signatures(synth.RD.repeats["X"], rms.event_unit(event), context, S)
    (signature,) = runs
    return signature


def _sites(records: list[Any], s: HybridSettings = S) -> list[dict[str, Any]]:
    members = trm._event_peak(records)
    sample = random.Random(s.seed).sample(members, min(len(members), s.phase_run_minor_max_reads))
    clean: list[dict[tuple[str, int], int]] = []
    impure: list[dict[tuple[str, int], int]] = []
    _f, meta = features(
        synth.allele(rms.LONG), [m.seq for m in sample], s, clean=clean, impure=impure
    )
    strands = [m.strand for m in sample]
    return run_minor_sites(clean, strands, meta, s, KNOWN, impure=impure)


@pytest.mark.parametrize("event", IN_RUN_EVENTS)
def test_signatures_shorten_the_x_unit_c7_run(event: str) -> None:
    base_, parent, carrier = _signature(event)
    assert (base_, parent) == ("C", 7) and carrier < parent


@pytest.mark.parametrize("event", IN_RUN_EVENTS)
def test_impure_observation_is_the_longest_stretch_inside_the_bounding_bases(event: str) -> None:
    cons = synth.allele(rms.LONG)
    read = rms.carrier(rms.LONG, event=event)
    start = base._unit_run(cons, rms.EVENT_UNIT)
    clean: list[dict[tuple[str, int], int]] = []
    impure: list[dict[tuple[str, int], int]] = []
    (f,), _meta = features(cons, [read], S, clean=clean, impure=impure)
    assert ("run", start) not in clean[0]
    assert impure[0][("run", start)] == f[("run", start)] == _signature(event)[2]
    # Every other run of the read is clean, so it has no impure observation.
    assert set(impure[0]) == {("run", start)}


@pytest.mark.parametrize("event", IN_RUN_EVENTS)
@pytest.mark.parametrize("shape", sorted(rms.SHAPES))
def test_in_run_event_minority_is_found(shape: str, event: str) -> None:
    records = rms.minority(shape, trm.AF, S.phase_run_minor_max_reads, SEED, "two_peak", event)
    sites = _sites(records)
    assert sites and sites[0]["site"] == trm._event_site(), sites[:3]
    assert (sites[0]["major"], sites[0]["minor"]) == _signature(event)[1:]
    assert sites[0]["observation"] == "in_run"
    assert sites[0]["share_lower_bound"] >= S.phase_run_minor_min_share


@pytest.mark.parametrize("shape", sorted(rms.SHAPES))
def test_in_run_noise_on_wild_type_is_expected(shape: str) -> None:
    records = rms.in_run_noise(shape, NOISE_RATE, S.phase_run_minor_max_reads, SEED, "two_peak")
    assert [s for s in _sites(records) if s["observation"] == "in_run"] == []


def test_in_run_test_can_be_switched_off_and_needs_a_signature() -> None:
    records = rms.minority("hifi", trm.AF, S.phase_run_minor_max_reads, SEED, "two_peak", "insG")
    off = dataclasses.replace(S, phase_run_minor_in_run=False)
    assert [s for s in _sites(records, off) if s["observation"] == "in_run"] == []
    scope_off = dataclasses.replace(S, phase_run_minor_scope="off")
    assert _sites(records, scope_off) == []
    members = trm._event_peak(records)
    clean: list[dict[tuple[str, int], int]] = []
    impure: list[dict[tuple[str, int], int]] = []
    _f, meta = features(
        synth.allele(rms.LONG), [m.seq for m in members], S, clean=clean, impure=impure
    )
    strands = [m.strand for m in members]
    no_signatures = dataclasses.replace(KNOWN, runs=frozenset())
    sites = run_minor_sites(clean, strands, meta, S, no_signatures, impure=impure)
    assert [s for s in sites if s["observation"] == "in_run"] == []


def test_in_run_setting_is_validated() -> None:
    with pytest.raises(ValueError, match="phase_run_minor_in_run"):
        HybridSettings(phase_run_minor_in_run="yes")  # type: ignore[arg-type]
