"""A site carrying a known event's site-table signature is never dropped.

Found by an adversarial sweep: in a two-peak sample, a real
minority haplotype carrying dupC (C7 -> C8) at an allele fraction of 0.15, whose
carriers are all low-quality reads, plus insertion stutter at that run in other
low-quality reads, reaches a run-site minor AF of about 0.25. The 15j low-accuracy
test then explains the site by poor reads (its high-quality minor AF is about 0.08)
and drops it, so the sample became NEGATIVE; with the rule off it is INCONCLUSIVE.

The quality rules (15j two-peak site drop, 15k single-event alternative) keep any site
whose change is the site-table signature of a dictionary template. The signatures are
derived by running every template (insertions, deletions, delete-inserts) in each
allowed unit through the site table: dupC lengthens the X unit's C7 run, insG and
delinsAT shorten it (the inserted G splits the run), the deletions shorten C3/C4 runs,
dupA adds an insertion slot. Column sites are never signatures; insG_pos54 in unit J
changed no site at all until the site table recorded the insertion slots next to a run.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.config import apply_mutation
from muc_one_span.hybrid.known_events import is_known_event_site, known_event_sites
from muc_one_span.hybrid.phase_quality import quality_sites
from muc_one_span.hybrid.phase_sites import candidates, features
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_quality_groups as groups
from tests.unit.hybrid import test_quality_single_event as single
from tests.unit.hybrid import test_quality_sites as quality
from tests.unit.hybrid import test_single_event as base

# The heavy synthetic sweeps below carry the safety_sweep marker (their own CI job and
# make test-unit, never skipped); the cheap tests run in the core suite.

S = quality.S  # the quality rules are opt-in since switched on explicitly
OPT_IN = quality.OPT_IN
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
FOUND = ("ont", 0.15, 1)  # (profile, minority AF, seed) of the shape found in the sweep
FLOOR_SHAPE = (0.20, 2)  # (minority AF, seed) of the 15j/15l shape below
SWEEP_AFS = (0.15, 0.20, 0.25)
SWEEP_SEEDS = (0, 1, 2, 3)


def _off() -> RuntimeSettings:
    return dataclasses.replace(
        DEFAULT_SETTINGS, hybrid=dataclasses.replace(S, phase_quality_alpha=0.0)
    )


KNOWN = known_event_sites(synth.RD, S)
CONTEXT = synth.RD.repeats[synth.RD.canonical_repeat]
# Template/unit pairs whose mutation changes no run or insertion slot of the site
# table. insG_pos54 in unit J was one until the slots next to a run were recorded (its G lands in the insertion slot
# right before a C run, which the site table did not record); none is left.
BLIND: tuple[tuple[str, str], ...] = ()
TEMPLATE_UNITS = [
    (name, unit)
    for name, template in synth.RD.mutations.items()
    for unit in template["allowed_repeats"]
    if unit in synth.RD.repeats
]


def _event_sites(name: str, unit: str) -> list[tuple[dict[str, object], dict[object, object]]]:
    """Each non-column site the template changes, as a within-peak minor (parent major)."""
    parent = synth.RD.repeats[unit]
    mutated = apply_mutation(parent, synth.RD.mutations[name]["changes"])
    cons = CONTEXT + parent + CONTEXT
    (ref, mut), meta = features(cons, [cons, CONTEXT + mutated + CONTEXT], S)
    return [
        ({"site": k, "major": ref.get(k), "minor": mut.get(k)}, meta)
        for k in set(ref) | set(mut)
        if k[0] != "col" and ref.get(k) != mut.get(k) and ref.get(k) is not None
    ]


@pytest.mark.parametrize(("name", "unit"), TEMPLATE_UNITS)
def test_every_template_signature_is_protected(name: str, unit: str) -> None:
    """Insertions, deletions and delete-inserts: every run or insertion-slot change a
    template makes in any allowed unit is protected (insG and delinsAT shorten the C7
    run; the deletions shorten C3/C4 runs)."""
    sites = _event_sites(name, unit)
    if (name, unit) in BLIND:
        assert not sites and (name, unit) in KNOWN.blind
        return
    assert sites, "the template must change a run or an insertion slot"
    for site, meta in sites:
        assert is_known_event_site(site, meta, KNOWN), (name, unit, site)  # type: ignore[arg-type]


def test_blind_templates_are_listed() -> None:
    """Every template changes a run or an insertion slot of the table."""
    assert KNOWN.blind == BLIND


@pytest.mark.parametrize(
    ("site", "meta", "protected"),
    [
        ({"site": ("run", 9), "major": 7, "minor": 8}, {("run", 9): ("C", 7)}, True),
        ({"site": ("run", 9), "major": 7, "minor": 6}, {("run", 9): ("C", 7)}, True),
        ({"site": ("run", 9), "major": 7, "minor": 4}, {("run", 9): ("C", 7)}, True),
        ({"site": ("run", 9), "major": 3, "minor": 2}, {("run", 9): ("C", 3)}, False),
        ({"site": ("ins", 9), "major": "", "minor": "A"}, {}, True),
        ({"site": ("ins", 9), "major": "", "minor": "TT"}, {}, False),
        ({"site": ("col", 9), "major": "G", "minor": "C"}, {}, False),
    ],
)
def test_only_template_signatures_are_protected(
    site: dict[str, object], meta: dict[object, object], protected: bool
) -> None:
    assert is_known_event_site(site, meta, KNOWN) is protected  # type: ignore[arg-type]


@pytest.mark.safety_sweep
def test_all_poor_dupc_minority_is_not_released_to_negative(tmp_path: Path) -> None:
    records = groups._minority(*FOUND, dupc=True)
    for sub in ("off", "on"):
        (tmp_path / sub).mkdir()
    _off_summary, off = base._run(tmp_path / "off", records, _off())
    assert off["state"] == "INCONCLUSIVE", "precondition: the site blocks NEGATIVE without the rule"
    summary, decision = base._run(tmp_path / "on", records, OPT_IN)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["quality_associated_sites"]
    assert not any(
        s["kind"] == "run" and s["minor"] > s["major"]
        for s in summary["hybrid"]["quality_associated_sites"]
    )


@pytest.mark.safety_sweep
def test_an_explained_substitution_is_still_dropped() -> None:
    """The guard is narrow: the 15j column artefact is still explained and dropped."""
    members = quality._long_peak(quality._sample(quality.SEEDS[0]))
    cons = synth.allele(quality.LONG)
    feats, meta = features(cons, [m.seq for m in members], S)
    sites = candidates(feats, [m.strand for m in members], meta, S)
    quals = [m.mean_q for m in members]
    kept, dropped = quality_sites(sites, feats, quals, S, meta=meta, insertions=KNOWN)
    assert dropped and not kept


# With every quality rule off (the defaults) the minority's dupC run and its
# substitutions stayed below the candidate tiers at 15% of one allele (the within-peak
# floor), and with the opt-in 15j rule FLOOR_SHAPE's only visible marker, one
# substitution site, was explained away: 16 strict xfails in 15k. The run-minority
# tier (hybrid.run_minor) now makes the dupC run visible and the known-event guard
# keeps it, so every case must pass, with and without the opt-in rule.
SWEEP = [(p, af, seed) for p in sorted(single.PROFILES) for af in SWEEP_AFS for seed in SWEEP_SEEDS]


@pytest.mark.safety_sweep
@pytest.mark.parametrize(("profile", "af", "seed"), SWEEP)
def test_all_poor_dupc_minority_is_never_negative_at_the_defaults(
    tmp_path: Path, profile: str, af: float, seed: int
) -> None:
    records = groups._minority(profile, af, seed, dupc=True)
    summary, decision = base._run(tmp_path, records)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["selection_detail"]


@pytest.mark.safety_sweep
@pytest.mark.parametrize(("profile", "af", "seed"), SWEEP)
def test_all_poor_dupc_minority_is_never_negative_opted_in(
    tmp_path: Path, profile: str, af: float, seed: int
) -> None:
    records = groups._minority(profile, af, seed, dupc=True)
    summary, decision = base._run(tmp_path, records, OPT_IN)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["quality_associated_sites"]


@pytest.mark.safety_sweep
@pytest.mark.parametrize("profile", sorted(single.PROFILES))
def test_floor_shape_is_inconclusive_at_the_defaults(tmp_path: Path, profile: str) -> None:
    """The 15l shape at the shipped defaults (quality rules off): the minority's
    substitution site stays and blocks a negative call."""
    summary, decision = base._run(tmp_path, groups._minority(profile, *FLOOR_SHAPE, dupc=True))
    assert decision["state"] == "INCONCLUSIVE", summary["hybrid"]["selection_detail"]


@pytest.mark.safety_sweep
def test_a_guard_kept_site_that_poor_reads_explain_is_reported() -> None:
    """``guarded_but_explained`` names a site kept only by the known-event
    guard although the 15j test explains it by poor reads (the guard widened to every
    site here, so the explained column artefact of the 15j shape is such a site)."""
    from unittest.mock import patch

    from muc_one_span.hybrid import phase_quality

    members = quality._long_peak(quality._sample(quality.SEEDS[0]))
    cons = synth.allele(quality.LONG)
    feats, meta = features(cons, [m.seq for m in members], S)
    sites = candidates(feats, [m.strand for m in members], meta, S)
    quals = [m.mean_q for m in members]
    args = (sites, feats, quals, S)
    assert phase_quality.guarded_but_explained(*args, meta=meta, insertions=KNOWN) == []
    with patch.object(phase_quality, "is_known_event_site", return_value=True):
        kept, dropped = quality_sites(*args, meta=meta, insertions=KNOWN)
        guarded = phase_quality.guarded_but_explained(*args, meta=meta, insertions=KNOWN)
    assert kept == sites and not dropped
    assert guarded == sites


def test_single_event_alternative_fails_closed_on_a_guard_only_site() -> None:
    """The opt-in single-event alternative is never offered when its one
    remaining event rests on a site kept only by the guard but explained by poor reads
    (for example a protected C7 -> C6 stutter artefact): the peak keeps its own result."""
    from types import SimpleNamespace
    from unittest.mock import patch

    from muc_one_span.hybrid import phase

    kept, dropped = [{"site": ("run", 1)}], [{"site": ("col", 2)}]
    offered = phase.PhaseResult([[]], phase.UNCONFIRMED_SINGLE_SITE)

    def table(guard_only: list[dict[str, object]]) -> SimpleNamespace:
        return SimpleNamespace(
            meta={},
            quality=lambda sites, settings: (kept, dropped),
            guard_only=lambda sites, settings: guard_only,
        )

    with (
        patch.object(phase, "events", side_effect=lambda sites, meta: [[s] for s in sites]),
        patch.object(phase, "_split", return_value=offered),
    ):
        both = kept + dropped
        alternative: Any = phase._single_event_alternative
        assert alternative(table([]), both, S) is offered
        assert alternative(table(kept), both, S) is None


def test_blind_template_pairs_are_warned_when_a_signature_consumer_is_on(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A custom dictionary's template/unit pair without a site signature can be
    protected neither by the quality rules' guard nor seen by the run-minority tier's
    signature tests (the known-event scope and the within-run test); a run with any of
    them on says so."""
    import logging

    from muc_one_span.hybrid.known_events import KnownEventSites, warn_blind_templates

    blind = KnownEventSites(blind=(("subst_x", "X"),))
    rules_off = dataclasses.replace(S, phase_quality_alpha=0.0)
    tier_off = dataclasses.replace(rules_off, phase_run_minor_scope="off")
    in_run_off = dataclasses.replace(rules_off, phase_run_minor_in_run=False)
    known_scope = dataclasses.replace(in_run_off, phase_run_minor_scope="known_events")
    with caplog.at_level(logging.WARNING, logger="muc_one_span.hybrid.known_events"):
        warn_blind_templates(blind, tier_off)
        warn_blind_templates(blind, in_run_off)  # scope "all": every run, no signatures
        assert not caplog.records  # no signature consumer is on
        warn_blind_templates(KNOWN, S)
        assert not caplog.records  # the bundled dictionary has no blind pair
        for settings in (S, rules_off, known_scope):
            warn_blind_templates(blind, settings)
    assert len(caplog.records) == 3
    assert all("subst_x in X" in r.getMessage() for r in caplog.records)
