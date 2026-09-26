"""Task 13b: read support must beat the read-derived alternative at the event site.

A consensus error (a POA/polish indel) must never be rated ``supported`` just because
reads match the consensus window better than the dictionary parent. Each event is
compared with the no-event allele (the consensus with the event unit reverted to its
parent) and with the most frequent read-derived alternative over a context that
extends ``event_context_units`` repeat units beyond the unit; a read that favours
neither side is not support, and an event whose estimated alternative share exceeds
``event_max_alternative_frac`` is not ``supported``.
"""

from __future__ import annotations

import json
import random
import re
from dataclasses import replace
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.clinical_gates import mutation_blockers, mutation_supported
from muc_one_span.hybrid.engine import annotate_read_support
from muc_one_span.hybrid.evidence import event_allele_fraction, event_read_support
from muc_one_span.hybrid.spans import Anchors, ReadRecord, categorize_reads
from muc_one_span.pipeline import execute_pipeline
from muc_one_span.report import compute_clinical_decision
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid.synth import LAYOUT

S = HybridSettings()
X = synth.RD.repeats["X"]
UNIT_I = synth.RD.repeats["I"]
LACKING = [0.4, 0.5, 0.6]  # share of reads that lack the consensus error
N_READS = 60
ERR = 0.02


def _stutter(seq: str, rng: random.Random, p_del: float, p_ins: float) -> str:
    """Per-read homopolymer stutter: each run >= hp_event_min_run gains/loses one base."""
    out = []
    for match in re.finditer(r"(.)\1*", seq):
        run = match.group(0)
        if len(run) >= S.hp_event_min_run:
            r = rng.random()
            run = run[:-1] if r < p_del else run + run[0] if r < p_del + p_ins else run
        out.append(run)
    return "".join(out)


def _oriented_mix(
    parts: list[tuple[str, int]], seed: int, *, p_del: float = 0.0, p_ins: float = 0.0
) -> list[tuple[str, str]]:
    """Spanning reads (oriented, with strand) from several source alleles."""
    rng = random.Random(seed)
    records: list[ReadRecord] = []
    for source, count in parts:
        for _ in range(count):
            template = _stutter(source, rng, p_del, p_ins)
            records += synth.reads(template, 1, err=ERR, seed=rng.randrange(1 << 30))
    cats = categorize_reads(records, Anchors.from_dictionary(synth.RD, S, LAYOUT), S)
    return [(sp.seq, sp.strand) for sp in cats.spanning]


def _annotated(cons: str, reads: list[tuple[str, str]], s: HybridSettings = S) -> dict[str, Any]:
    result = classify_sequence(cons, synth.RD)
    return annotate_read_support(
        "allele_1", cons, result, rd=synth.RD, members={"allele_1": reads}, settings=s
    )


def _split(n_lacking: float) -> tuple[int, int]:
    lacking = round(N_READS * n_lacking)
    return N_READS - lacking, lacking


# --- (a) a 1-base consensus indel that many reads lack ---------------------------------


@pytest.mark.parametrize("lacking", LACKING)
def test_mid_unit_consensus_insertion_is_not_supported(lacking: float) -> None:
    mid = len(X) // 2
    error_unit = X[:mid] + "T" + X[mid:]
    truth = synth.allele(["X"] * 8 + ["X"] + ["X"] * 8)
    cons = synth.allele(["X"] * 8 + [error_unit] + ["X"] * 8)
    with_error, without = _split(lacking)
    reads = _oriented_mix([(cons, with_error), (truth, without)], seed=31)
    result = _annotated(cons, reads)
    assert result["mutations_detected"], "fixture: the injected indel must be detected"
    for mutation in result["mutations_detected"]:
        support = mutation["read_support"]
        assert support["status"] != "supported", support
        assert not mutation_supported(mutation)


@pytest.mark.parametrize("lacking", LACKING)
def test_boundary_displaced_consensus_error_is_not_supported(lacking: float) -> None:
    """The M3_hifi shape: an extra C in unit I's C5 run is reported one unit later."""
    run5 = "G" + "C" * 5 + "A"
    assert run5 in UNIT_I  # fixture: dictionary unit I ends GCG C5 A
    error_i = UNIT_I.replace(run5, "G" + "C" * 6 + "A")
    truth = synth.allele(["X"] * 8 + ["I", "X"] + ["X"] * 8)
    cons = synth.allele(["X"] * 8 + [error_i, "X"] + ["X"] * 8)
    with_error, without = _split(lacking)
    reads = _oriented_mix([(cons, with_error), (truth, without)], seed=32)
    result = _annotated(cons, reads)
    assert result["mutations_detected"], "fixture: the consensus error must be detected"
    for mutation in result["mutations_detected"]:
        assert mutation["read_support"]["status"] != "supported", mutation


# --- (b) an extra C at the dupC position on a wild-type sample -------------------------


@pytest.mark.parametrize("seed", range(40, 45))
@pytest.mark.parametrize("lacking", LACKING)
def test_consensus_dupc_on_a_wild_type_mixture_is_never_pathogenic(
    lacking: float, seed: int
) -> None:
    wild = synth.allele(["X"] * 10 + ["X"] + ["X"] * 10)
    cons = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    with_error, without = _split(lacking)
    reads = _oriented_mix([(cons, with_error), (wild, without)], seed=seed, p_del=0.1, p_ins=0.05)
    result = _annotated(cons, reads)
    dupc = [m for m in result["mutations_detected"] if m.get("mutation_name") == "dupC"]
    assert len(dupc) == 1
    assert dupc[0]["read_support"]["kind"] == "homopolymer"
    assert dupc[0]["read_support"]["status"] != "supported", dupc[0]["read_support"]
    assert any(b.startswith("read-level support") for b in mutation_blockers(dupc[0]))


# --- (c) true events at realistic stutter stay supported -------------------------------


@pytest.mark.parametrize(("p_del", "p_ins"), [(0.1, 0.05), (0.25, 0.08)])
def test_true_dupc_with_realistic_stutter_stays_supported(p_del: float, p_ins: float) -> None:
    cons = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    reads = _oriented_mix([(cons, N_READS)], seed=34, p_del=p_del, p_ins=p_ins)
    result = _annotated(cons, reads)
    dupc = next(m for m in result["mutations_detected"] if m.get("mutation_name") == "dupC")
    support = dupc["read_support"]
    assert support["status"] == "supported", support
    assert support["alternative_frac"] <= S.event_max_alternative_frac


def test_true_competition_event_stays_supported() -> None:
    dupa = next(
        seq for seq, (p, n) in synth.RD.mutated_sequences.items() if p == "X" and n == "dupA"
    )
    cons = synth.allele(["X"] * 10 + [dupa] + ["X"] * 10)
    reads = _oriented_mix([(cons, N_READS)], seed=35, p_del=0.1, p_ins=0.05)
    result = _annotated(cons, reads)
    event = next(m for m in result["mutations_detected"] if m.get("mutation_name") == "dupA")
    assert event["read_support"]["kind"] == "competition"
    assert event["read_support"]["status"] == "supported", event["read_support"]


# --- (d) heterozygous events at 3:1 imbalance stay supported ---------------------------


def _fastq(path: Path, reads: list[ReadRecord]) -> Path:
    path.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in reads))
    return path


@pytest.mark.parametrize(("n_normal", "n_carrier"), [(150, 50), (50, 150)])
def test_heterozygous_dupc_at_three_to_one_imbalance_stays_pathogenic(
    tmp_path: Path, n_normal: int, n_carrier: int
) -> None:
    normal = synth.allele(["X"] * 25)
    carrier = synth.allele(["X"] * 14 + [synth.dupc()] + ["X"] * 30)
    reads = synth.reads(normal, n_normal, err=ERR, seed=36) + synth.reads(
        carrier, n_carrier, err=ERR, seed=37
    )
    out = tmp_path / "out"
    with patch("muc_one_span.tools.check_tools"):
        execute_pipeline(
            str(_fastq(tmp_path / "in.fastq", reads)),
            str(out),
            None,
            "",
            1,
            10,
            5.0,
            False,
            "ont",
            None,
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
        )
    summary = json.loads((out / "summary.json").read_text())
    mutations = [m for c in summary["classifications"].values() for m in c["mutations"]]
    dupc = [m for m in mutations if m["mutation_name"] == "dupC"]
    assert len(dupc) == 1 and dupc[0]["read_support"]["status"] == "supported", dupc
    assert compute_clinical_decision(summary)["state"] == "PATHOGENIC"


# --- the new tunables control behaviour ------------------------------------------------


def test_max_alternative_frac_setting_controls_the_verdict() -> None:
    mid = len(X) // 2
    cons = synth.allele(["X"] * 8 + [X[:mid] + "T" + X[mid:]] + ["X"] * 8)
    truth = synth.allele(["X"] * 17)
    with_error, without = _split(min(LACKING))
    reads = _oriented_mix([(cons, with_error), (truth, without)], seed=38)
    strict = _annotated(cons, reads)["mutations_detected"][0]["read_support"]
    assert strict["status"] != "supported"
    lax = replace(S, event_max_alternative_frac=1.0)
    loose = _annotated(cons, reads, lax)["mutations_detected"][0]["read_support"]
    assert loose["status"] == "supported", loose
    assert loose["alternative_frac"] == strict["alternative_frac"] > S.event_max_alternative_frac


def test_context_units_setting_controls_the_compared_window() -> None:
    """Without context the displaced M3 error is re-rated supported (the old defect)."""
    run5 = "G" + "C" * 5 + "A"
    error_i = UNIT_I.replace(run5, "G" + "C" * 6 + "A")
    truth = synth.allele(["X"] * 8 + ["I", "X"] + ["X"] * 8)
    cons = synth.allele(["X"] * 8 + [error_i, "X"] + ["X"] * 8)
    with_error, without = _split(max(LACKING))
    reads = _oriented_mix([(cons, with_error), (truth, without)], seed=39)
    classification = classify_sequence(cons, synth.RD)
    default = event_read_support(cons, classification, reads, synth.RD, S)
    assert default and all(v["status"] != "supported" for v in default.values())
    unit_only = replace(S, event_context_units=0.0)
    narrow = event_read_support(cons, classification, reads, synth.RD, unit_only)
    assert any(v["status"] == "supported" for v in narrow.values()), narrow


# --- mixture fraction ------------------------------------------------------------------


def test_event_allele_fraction_recovers_the_mixture_weight() -> None:
    p_event, p_none = (0.8, 0.1), (0.1, 0.8)  # (p under event, p under no-event)
    assert event_allele_fraction([p_event] * 10) == 1.0
    assert event_allele_fraction([p_none] * 10) == 0.0
    half = event_allele_fraction([p_event] * 10 + [p_none] * 10)
    assert half == pytest.approx(0.5, abs=1e-9)
    assert event_allele_fraction([]) == 0.0  # no information fails closed
    assert event_allele_fraction([(0.3, 0.3)] * 5) == 0.0
