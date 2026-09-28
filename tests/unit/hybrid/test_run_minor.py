"""Task 15l: a within-peak run-length minority above the expected stutter blocks NEGATIVE.

The shape (found in the Task 15j fix round): two length peaks, and inside one allele's
peak a real minority haplotype (mosaicism, a third haplotype, a chimera) carrying dupC
(an X unit's C7 run read as C8) at 15-30% of that allele's reads, with realistic run
stutter. Its C8 share stays below every earlier run floor (``phase_run_bg_multiplier``
and the single-peak-only 15g ``phase_run_safety_multiplier``, both multiples of the
peers' raw share), so the sample was NEGATIVE.

The run-minority tier (``run_minor``) explains each run's clean observations as a
mixture of its modal length and another length, each convolved with the Task 15f
stutter profile of its own length from the peak's peer runs, and blocks a negative call
when the one-sided lower confidence bound of the minority share reaches
``phase_run_minor_min_share``. It never splits a peak or creates an event.

Also here (Task 15l addendum): an insertion slot next to a run is recorded (insG_pos54
in unit J inserts a G right before a C run and was invisible to the site table).
"""

from __future__ import annotations

import dataclasses
import random
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.config import _apply_mutation
from muc_one_span.hybrid.engine import RUN_MINOR_BASIS, reconstruct_alleles
from muc_one_span.hybrid.known_events import KnownEventSites, known_event_sites
from muc_one_span.hybrid.lengths import fit_length_model
from muc_one_span.hybrid.phase_sites import features
from muc_one_span.hybrid.run_minor import run_minor_sites
from muc_one_span.hybrid.spans import ReadRecord, SpanRead, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import run_minor_synth as rms
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_single_event as base

S = DEFAULT_SETTINGS.hybrid
KNOWN = known_event_sites(synth.RD, S)
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
# Site-level tests: the tier's full read sample per allele peak and a minority at the
# upper end of the documented floors (docs/reference/limitations.md). End-to-end
# tests: the phase-sample depth per allele and a minority well above every floor.
SITE_DEPTH = S.phase_run_minor_max_reads
AF = 0.30
DEPTH = 2 * S.phase_max_site_reads
E2E_AF = 0.40
SEEDS = (0, 1)


def _event_peak(records: list[ReadRecord]) -> list[SpanRead]:
    """Members of the event allele's (longest) length peak."""
    spans = categorize_reads(records, base.ANCH, S).spanning
    model = fit_length_model(spans, S, base.ANCH)
    return max(model.peaks, key=lambda p: p.center_bp).members


def _table(
    members: list[SpanRead], s: HybridSettings
) -> tuple[list[dict[tuple[str, int], int]], list[str], dict[tuple[str, int], tuple[str, int]]]:
    """(clean run observations, strands, run meta) of the phase sample of ``members``."""
    sample = random.Random(s.seed).sample(members, min(len(members), s.phase_max_site_reads))
    clean: list[dict[tuple[str, int], int]] = []
    _feats, meta = features(synth.allele(rms.LONG), [m.seq for m in sample], s, clean=clean)
    return clean, [m.strand for m in sample], meta


def _sites(members: list[SpanRead], s: HybridSettings = S) -> list[dict[str, Any]]:
    clean, strands, meta = _table(members, s)
    return run_minor_sites(clean, strands, meta, s, KNOWN)


def _event_site() -> tuple[str, int]:
    return ("run", base._unit_run(synth.allele(rms.LONG), rms.EVENT_UNIT))


# --- the site-level test ---------------------------------------------------------------


@pytest.mark.parametrize("shape", sorted(rms.SHAPES))
def test_minority_run_is_found_above_its_expected_stutter(shape: str) -> None:
    members = _event_peak(rms.minority(shape, AF, SITE_DEPTH, SEEDS[0], "two_peak"))
    sites = _sites(members)
    assert sites and sites[0]["site"] == _event_site(), sites[:3]
    assert (sites[0]["major"], sites[0]["minor"]) == (7, 8)
    assert sites[0]["share_lower_bound"] >= S.phase_run_minor_min_share


@pytest.mark.parametrize("shape", sorted(rms.SHAPES))
def test_wild_type_stutter_is_expected(shape: str) -> None:
    """Every run stutters by its length and strand, poor reads twice as much: no site."""
    members = _event_peak(rms.minority(shape, 0.0, SITE_DEPTH, SEEDS[0], "two_peak"))
    assert _sites(members) == []


def test_scope_off_and_known_events() -> None:
    members = _event_peak(rms.minority("ont", AF, SITE_DEPTH, SEEDS[0], "two_peak"))
    off = dataclasses.replace(S, phase_run_minor_scope="off")
    assert _sites(members, off) == []
    known = dataclasses.replace(S, phase_run_minor_scope="known_events")
    assert [s["site"] for s in _sites(members, known)][:1] == [_event_site()]


def test_known_events_scope_tests_only_template_signatures() -> None:
    """With no dictionary signature (C7 -> C8 is dupC's) the known-event scope is silent."""
    members = _event_peak(rms.minority("ont", AF, SITE_DEPTH, SEEDS[0], "two_peak"))
    known = dataclasses.replace(S, phase_run_minor_scope="known_events")
    clean, strands, meta = _table(members, known)
    assert run_minor_sites(clean, strands, meta, known, KnownEventSites()) == []


def test_clean_observation_skips_a_changed_bounding_base() -> None:
    """``CCCCCC[A]A`` read as ``CCCCCCCA`` (a neighbouring substitution, the simulated
    HiFi systematic error): the site table sees a C7 run, the clean observation none."""
    cons = synth.allele(rms.ARTEFACT_WT)
    start = base._unit_run(cons, 0)
    end = start + 6
    read = cons[:end] + "C" + cons[end + 1 :]
    clean: list[dict[tuple[str, int], int]] = []
    (f,), _meta = features(cons, [read], S, clean=clean)
    assert f[("run", start)] == 7
    assert ("run", start) not in clean[0]
    inserted = cons[:end] + "C" + cons[end:]
    clean = []
    features(cons, [inserted], S, clean=clean)
    assert clean[0][("run", start)] == 7


# --- end to end ------------------------------------------------------------------------


@pytest.mark.parametrize("seed", SEEDS)
def test_two_peak_run_minority_is_inconclusive_and_located(tmp_path: Path, seed: int) -> None:
    summary, decision = base._run(tmp_path, rms.minority("hifi", E2E_AF, DEPTH, seed, "two_peak"))
    assert decision["state"] == "INCONCLUSIVE", decision["details"]
    assert RUN_MINOR_BASIS in summary["hybrid"]["split_bases"]
    assert summary["hybrid"]["selection_status"] == "unresolved_run_minor"
    text = f"unresolved heterozygous site at repeat {rms.EVENT_REPEAT} (run 7>8"
    assert any(text in d for d in decision["details"]), decision["details"]


def test_run_minority_tier_never_splits_or_creates_an_event(tmp_path: Path) -> None:
    records = rms.minority("ont", E2E_AF, DEPTH, SEEDS[0], "two_peak")
    result = reconstruct_alleles(
        base._fastq(tmp_path / "in.fastq", records), tmp_path, synth.RD, DEFAULT_SETTINGS
    )
    assert RUN_MINOR_BASIS in result.block["split_bases"]
    assert "allele_2" in result.alleles and "allele_3" not in result.alleles
    for name in ("allele_1", "allele_2"):
        assert synth.dupc() not in result.consensus_paths[name].read_text()


@pytest.mark.parametrize("seed", SEEDS)
def test_single_peak_run_minority_is_never_negative(tmp_path: Path, seed: int) -> None:
    summary, decision = base._run(tmp_path, rms.minority("hifi", E2E_AF, DEPTH, seed, "one_peak"))
    assert decision["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


# --- insertion slots next to a run (insG_pos54 in unit J) ------------------------------

J_EVENT_UNIT = 10


def _j_alleles() -> tuple[str, str]:
    wild = ["X"] * J_EVENT_UNIT + ["J"] + ["X"] * (len(rms.LONG) - J_EVENT_UNIT - 1)
    mutated = _apply_mutation(synth.RD.repeats["J"], synth.RD.mutations["insG_pos54"]["changes"])
    carrier = [*wild[:J_EVENT_UNIT], mutated, *wild[J_EVENT_UNIT + 1 :]]
    return synth.allele(wild), synth.allele(carrier)


def test_insertion_before_a_run_is_a_site() -> None:
    wild, carrier = _j_alleles()
    feats, _meta = features(wild, [wild, carrier], S)
    changed = {k for k in feats[0] if feats[0][k] != feats[1].get(k)}
    assert any(k[0] == "ins" and feats[1][k] == "G" for k in changed), changed


def test_run_lengthening_is_not_an_insertion_allele() -> None:
    """A base of the run itself inserted next to it lengthens the run (its run site)."""
    cons = synth.allele(rms.LONG)
    start = base._unit_run(cons, rms.EVENT_UNIT)
    longer = cons[:start] + "C" + cons[start:]
    feats, _meta = features(cons, [cons, longer], S)
    assert {k for k in feats[0] if feats[0][k] != feats[1].get(k)} == {("run", start)}


@pytest.mark.parametrize("seed", SEEDS)
def test_equal_length_ins_g_pos54_heterozygote_is_never_negative(tmp_path: Path, seed: int) -> None:
    """Both alleles have the same unit count; only the J unit's inserted G differs."""
    wild, carrier = _j_alleles()
    rng = random.Random(seed)
    half = S.phase_max_site_reads
    records = [
        rms.read(carrier if i < half else wild, "hifi", rng, f"j{i}") for i in range(2 * half)
    ]
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


# --- settings ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("phase_run_minor_scope", "some"),
        ("phase_run_minor_alpha", 0.0),
        ("phase_run_minor_alpha", 1.0),
        ("phase_run_minor_min_share", 0.0),
        ("phase_run_minor_min_share", 1.0),
    ],
)
def test_run_minor_settings_are_validated(key: str, value: object) -> None:
    with pytest.raises(ValueError, match=key):
        HybridSettings(**{key: value})  # type: ignore[arg-type]
