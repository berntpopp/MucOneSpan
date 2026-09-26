"""Task 15f: the homopolymer mixture fit uses a run-length-aware stutter background.

Deletion stutter grows with run length. A true dupC run (C8) therefore shows more C7
reads than the C7 background runs show C6 reads; the pre-15f model (the no-event
profile shifted by one base) assigned that excess to the no-event allele, so a true
dupC looked like a C7/C8 mixture. Each allele of the mixture is now convolved with
the stutter profile of its own run length: measured at peer runs of that length when
there are enough, otherwise extrapolated from the two nearest measured lengths.
"""

from __future__ import annotations

import dataclasses
import random
import re
from collections.abc import Callable
from typing import Any

import pytest

from muc_one_span.clinical_gates import mutation_blockers
from muc_one_span.hybrid.polish import _runs
from muc_one_span.hybrid.spans import Anchors, ReadRecord, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_evidence_alternative as alt
from tests.unit.hybrid.synth import LAYOUT

S = HybridSettings()
N_READS = 200
SEEDS = range(50, 53)
# The dupC run: the dictionary X unit's final C run plus one base.
EVENT_LEN = max(e - s for s, e, b in _runs(synth.dupc(), 2) if b == "C")
# Realistic length-dependent deletion stutter, from the development panels: about 31%
# of reads show C7 at a true C8 run (simulated HiFi D1/HD4), and the deletion share
# grows about 2.5-fold per base between adjacent measured run lengths (C6 -> C7 in
# simulated HiFi and in real ONT "+" reads). Insertion stutter is length-independent.
DEL_AT_EVENT = 0.31
GROWTH = 2.5
P_INS = 0.05
# The D1 shape: C7 runs show only 8.2% deletion stutter, so the step from C7 to the
# true C8 run (~3.8-fold) is steeper than any trend below it.
D1_DEL_BELOW_EVENT = 0.082


def _log_linear(length: int) -> float:
    return float(DEL_AT_EVENT * GROWTH ** (length - EVENT_LEN))


def _d1_step(length: int) -> float:
    if length >= EVENT_LEN:
        return DEL_AT_EVENT
    return float(D1_DEL_BELOW_EVENT * GROWTH ** (length - EVENT_LEN + 1))


def _stutter(seq: str, rng: random.Random, p_del: Callable[[int], float]) -> str:
    """Each run of two or more bases loses one base at its length's deletion share."""
    out = []
    for match in re.finditer(r"(.)\1*", seq):
        run = match.group(0)
        if len(run) > 1:
            r = rng.random()
            d = p_del(len(run))
            run = run[:-1] if r < d else run + run[0] if r < d + P_INS else run
        out.append(run)
    return "".join(out)


def _reads(
    parts: list[tuple[str, int]], seed: int, p_del: Callable[[int], float] = _log_linear
) -> list[tuple[str, str]]:
    rng = random.Random(seed)
    records: list[ReadRecord] = []
    for source, count in parts:
        for _ in range(count):
            template = _stutter(source, rng, p_del)
            records += synth.reads(template, 1, err=alt.ERR, seed=rng.randrange(1 << 30))
    cats = categorize_reads(records, Anchors.from_dictionary(synth.RD, S, LAYOUT), S)
    return [(sp.seq, sp.strand) for sp in cats.spanning]


WILD = synth.allele(["X"] * 10 + ["X"] + ["X"] * 10)
DUPC = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)


def _dupc_support(cons: str, reads: list[tuple[str, str]], s: HybridSettings = S) -> Any:
    result = alt._annotated(cons, reads, s)
    dupc = [m for m in result["mutations_detected"] if m.get("mutation_name") == "dupC"]
    return result, dupc


@pytest.mark.parametrize("seed", SEEDS)
def test_true_dupc_with_length_dependent_stutter_is_supported(seed: int) -> None:
    _result, dupc = _dupc_support(DUPC, _reads([(DUPC, N_READS)], seed))
    support = dupc[0]["read_support"]
    assert support["status"] == "supported", support
    assert support["alternative_frac"] <= S.event_max_alternative_frac / 2, support


@pytest.mark.parametrize("seed", SEEDS)
def test_true_dupc_with_a_d1_shaped_stutter_step_is_supported(seed: int) -> None:
    """The D1/HD4 shape (pre-15f: alternative share 0.23-0.29 here, 0.263/0.284 there)."""
    reads = _reads([(DUPC, N_READS)], seed, _d1_step)
    _result, dupc = _dupc_support(DUPC, reads)
    support = dupc[0]["read_support"]
    assert support["status"] == "supported", support
    legacy = dataclasses.replace(S, hp_stutter_model="shift")
    _result, old = _dupc_support(DUPC, reads, legacy)
    assert support["alternative_frac"] < old[0]["read_support"]["alternative_frac"]


@pytest.mark.parametrize("seed", SEEDS)
def test_wild_type_at_the_same_stutter_has_no_event(seed: int) -> None:
    reads = _reads([(WILD, N_READS)], seed)
    result, dupc = _dupc_support(WILD, reads)
    assert not dupc, result["mutations_detected"]
    # A spurious C8 consensus over the same wild-type reads is never supported.
    _result, spurious = _dupc_support(DUPC, reads)
    assert spurious[0]["read_support"]["status"] != "supported", spurious[0]["read_support"]
    assert any(b.startswith("read-level support") for b in mutation_blockers(spurious[0]))


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("wild_share", [0.4, 0.5])
def test_wild_type_dupc_mixture_is_discordant(wild_share: float, seed: int) -> None:
    wild = round(N_READS * wild_share)
    reads = _reads([(DUPC, N_READS - wild), (WILD, wild)], seed)
    _result, dupc = _dupc_support(DUPC, reads)
    support = dupc[0]["read_support"]
    # The mixture fit identifies the wild-type share (the discordant criterion); a
    # mixture whose log-likelihood ratio already fails is not_supported first.
    assert support["status"] in ("discordant", "not_supported"), support
    assert support["alternative_frac"] > S.event_max_alternative_frac, support


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("hp_stutter_model", "linear"),
        ("hp_stutter_min_class_runs", 0),
        ("hp_stutter_min_class_reads", 0),
        ("hp_stutter_max_growth", 0.5),
    ],
)
def test_stutter_settings_are_validated(field: str, value: object) -> None:
    with pytest.raises(ValueError, match=field):
        dataclasses.replace(S, **{field: value})


@pytest.mark.parametrize("field", ["hp_stutter_min_class_runs", "hp_stutter_min_class_reads"])
def test_sparse_length_classes_fall_back_to_the_shift_model(field: str) -> None:
    """With no length class measurable the fit is exactly the pre-15f shifted model."""
    reads = _reads([(DUPC, N_READS)], SEEDS[0])
    _r, legacy = _dupc_support(DUPC, reads, dataclasses.replace(S, hp_stutter_model="shift"))
    # More peer runs (or reads) than the allele can hold: no class is measured.
    sparse = dataclasses.replace(S, **{field: len(DUPC) * len(reads)})
    _r, dupc = _dupc_support(DUPC, reads, sparse)
    assert dupc[0]["read_support"] == legacy[0]["read_support"]
