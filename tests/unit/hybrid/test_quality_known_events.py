"""Task 15k: a site whose minor allele is a known insertion event is never dropped.

Found while building the Task 15k adversarial sweep: in a two-peak sample, a real
minority haplotype carrying dupC (C7 -> C8) at an allele fraction of 0.15, whose
carriers are all low-quality reads, plus insertion stutter at that run in other
low-quality reads, reaches a run-site minor AF of about 0.25. The 15j low-accuracy
test then explains the site by poor reads (its high-quality minor AF is about 0.08)
and drops it, so the sample became NEGATIVE; with the rule off it is INCONCLUSIVE.

The quality rules (15j two-peak site drop, 15k single-event alternative) now keep any
site whose minor allele adds exactly the inserted sequence of a dictionary insertion
template (dupC, dupA, insG, ...), taken from the bundled repeat dictionary. Other sites
(substitutions, deletions, insertions no template makes) are unaffected.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from muc_one_span.hybrid.phase_quality import inserts_known_event, known_insertions, quality_sites
from muc_one_span.hybrid.phase_sites import candidates, features
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_quality_groups as groups
from tests.unit.hybrid import test_quality_single_event as single
from tests.unit.hybrid import test_quality_sites as quality
from tests.unit.hybrid import test_single_event as base

S = DEFAULT_SETTINGS.hybrid
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
FOUND = ("ont", 0.15, 1)  # (profile, minority AF, seed) of the shape found in the sweep
SWEEP_AFS = (0.15, 0.20, 0.25)
SWEEP_SEEDS = (0, 1, 2, 3)


def _off() -> RuntimeSettings:
    return dataclasses.replace(
        DEFAULT_SETTINGS, hybrid=dataclasses.replace(S, phase_quality_alpha=0.0)
    )


def test_known_insertions_come_from_the_dictionary() -> None:
    known = known_insertions(synth.RD)
    templates = [
        str(c["sequence"]).upper()
        for t in synth.RD.mutations.values()
        for c in t["changes"]
        if c["type"] == "insert"
    ]
    assert known.inserted == frozenset(templates)
    # dupC lengthens the X unit's C7 run by one C; insC_pos23 the A unit's C4 run.
    assert ("C", 7, 1) in known.runs and ("C", 4, 1) in known.runs
    assert ("C", 3, 1) not in known.runs  # no template inserts into a C3 run


@pytest.mark.parametrize(
    ("site", "meta", "protected"),
    [
        ({"site": ("run", 9), "major": 7, "minor": 8}, {("run", 9): ("C", 7)}, True),
        ({"site": ("run", 9), "major": 3, "minor": 4}, {("run", 9): ("C", 3)}, False),
        ({"site": ("run", 9), "major": 7, "minor": 6}, {("run", 9): ("C", 7)}, False),
        ({"site": ("ins", 9), "major": "", "minor": "G"}, {}, True),
        ({"site": ("ins", 9), "major": "", "minor": "T"}, {}, False),
        ({"site": ("col", 9), "major": "G", "minor": "C"}, {}, False),
    ],
)
def test_only_template_insertions_are_protected(
    site: dict[str, object], meta: dict[object, object], protected: bool
) -> None:
    assert inserts_known_event(site, meta, known_insertions(synth.RD)) is protected  # type: ignore[arg-type]


def test_all_poor_dupc_minority_is_not_released_to_negative(tmp_path: Path) -> None:
    records = groups._minority(*FOUND, dupc=True)
    for sub in ("off", "on"):
        (tmp_path / sub).mkdir()
    _off_summary, off = base._run(tmp_path / "off", records, _off())
    assert off["state"] == "INCONCLUSIVE", "precondition: the site blocks NEGATIVE without the rule"
    summary, decision = base._run(tmp_path / "on", records)
    assert decision["state"] != NEGATIVE, summary["hybrid"]["quality_associated_sites"]
    assert not any(
        s["kind"] == "run" and s["minor"] > s["major"]
        for s in summary["hybrid"]["quality_associated_sites"]
    )


def test_an_explained_substitution_is_still_dropped() -> None:
    """The guard is narrow: the 15j column artefact is still explained and dropped."""
    members = quality._long_peak(quality._sample(quality.SEEDS[0]))
    cons = synth.allele(quality.LONG)
    feats, meta = features(cons, [m.seq for m in members], S)
    sites = candidates(feats, [m.strand for m in members], meta, S)
    quals = [m.mean_q for m in members]
    kept, dropped = quality_sites(
        sites, feats, quals, S, meta=meta, insertions=known_insertions(synth.RD)
    )
    assert dropped and not kept


# Pre-existing (Task 15j, not a 15k rule): here the minority's dupC run never becomes
# a candidate or safety-tier site (the within-peak detection floor under repair in
# Task 15l), so its only visible marker is one substitution site. All its carriers are
# low quality, so the 15j two-peak test explains that site away and the sample is
# NEGATIVE; with every quality rule off it is INCONCLUSIVE. Strict: 15l must revisit.
FLOOR_15L = pytest.mark.xfail(
    strict=True, reason="15j drop of the only marker of an invisible dupC minority (15l)"
)
SWEEP = [
    pytest.param(profile, af, seed, marks=FLOOR_15L)
    if (af, seed) == (0.20, 2)
    else (profile, af, seed)
    for profile in sorted(single.PROFILES)
    for af in SWEEP_AFS
    for seed in SWEEP_SEEDS
]


@pytest.mark.parametrize(("profile", "af", "seed"), SWEEP)
def test_quality_rules_never_release_an_all_poor_dupc_minority(
    tmp_path: Path, profile: str, af: float, seed: int
) -> None:
    """Never NEGATIVE unless it is NEGATIVE with every quality rule off as well (the
    pre-existing within-peak detection floor below the candidate tiers)."""
    records = groups._minority(profile, af, seed, dupc=True)
    for sub in ("off", "on"):
        (tmp_path / sub).mkdir()
    summary, decision = base._run(tmp_path / "on", records)
    if decision["state"] == NEGATIVE:
        _s, off = base._run(tmp_path / "off", records, _off())
        assert off["state"] == NEGATIVE, summary["hybrid"]["quality_associated_sites"]
