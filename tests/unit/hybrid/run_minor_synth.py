"""Synthetic within-peak run-minority samples (Task 15l tests and sweeps).

A sample holds one or two alleles; inside one allele's length peak a minority of
reads carries a run-length event (dupC: an X unit's C7 run read as C8) at allele
fraction ``af`` of that allele's reads. Every homopolymer run of every read stutters
by strand and run length (``SHAPES``: HiFi-like length-dependent stutter, ONT-like
strand-asymmetric stutter, and ONT-like saturating "+" strand stutter). Low-quality
reads (``POOR_FRAC``, qualities from ``QUALITIES``) stutter ``POOR_STUTTER_FACTOR``
times as often (quality-correlated stutter) and carry more random error.

``artefact`` samples are wild-type samples with a site-specific +1 run artefact at
one run (the simulated HiFi shape of a homozygous normal: one C unit's C6 run read as
C7 in about a third of the reads). A C insertion there makes the C unit read exactly
like an X unit carrying dupA, so the artefact reads are the carrier reads of a dupA
minority at the same share.
"""

from __future__ import annotations

import random
import re
from collections.abc import Callable

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.spans import PHRED_OFFSET, ReadRecord
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_run_site_safety as rss
from tests.unit.hybrid import test_single_event as base

Profile = Callable[[int], tuple[float, float]]
SHAPES: dict[str, dict[str, Profile]] = {
    "hifi": rss.SHAPES["hifi"],
    "ont": rss.SHAPES["ont_asymmetric"],
    "ont_saturating": rss.SHAPES["ont_saturating_plus"],
}
# Per-read mean Phred ranges by read class (ONT-like: disjoint; HiFi-like: overlapping).
QUALITIES = {
    "hifi": {"good": (30, 40), "poor": (22, 33)},
    "ont": {"good": (18, 26), "poor": (8, 17)},
    "ont_saturating": {"good": (18, 26), "poor": (8, 17)},
}
POOR_FRAC = 0.4
POOR_STUTTER_FACTOR = 2
POOR_ERR_FACTOR = 2
SHORT = ["X"] * 20
LONG = ["X"] * 30
EVENT_UNIT = 12  # 0-based inner unit of LONG whose C7 run carries the minority dupC
# 1-based repeat number of the event unit (the located reason names it).
EVENT_REPEAT = len(synth.PRE) + EVENT_UNIT + 1
LAYOUTS = ("two_peak", "one_peak")
# The artefact sample: a wild-type allele with one C unit (C6AA) among X units.
ARTEFACT_WT = ["C"] + ["X"] * 29


def _stutter(seq: str, rng: random.Random, profile: Profile, factor: float) -> str:
    """Each run of >= phase_run_min_len bases loses or gains one base by its length."""
    out = []
    for match in re.finditer(r"(.)\1*", seq):
        run = match.group(0)
        if len(run) >= rss.S.phase_run_min_len:
            p_del, p_ins = (min(p * factor, 1 / 2) for p in profile(len(run)))
            r = rng.random()
            run = run[:-1] if r < p_del else run + run[0] if r < p_del + p_ins else run
        out.append(run)
    return "".join(out)


def read(template: str, shape: str, rng: random.Random, name: str) -> ReadRecord:
    """One read of ``template``: random strand and quality class, stutter, error."""
    strand = "+" if rng.random() < 1 / 2 else "-"
    poor = rng.random() < POOR_FRAC
    factor = POOR_STUTTER_FACTOR if poor else 1
    stuttered = _stutter(template, rng, SHAPES[shape][strand], factor)
    err = base.ERR * (POOR_ERR_FACTOR if poor else 1)
    seq = synth.reads(stuttered, 1, err=err, seed=rng.randrange(1 << 30), strand_mix=False)[0].seq
    seq = seq if strand == "+" else rc(seq)
    q = QUALITIES[shape]["poor" if poor else "good"]
    return ReadRecord(name, seq, chr(rng.randint(*q) + PHRED_OFFSET) * len(seq))


def carrier(inner: list[str], unit: int = EVENT_UNIT) -> str:
    """``inner`` with dupC in ``unit``."""
    return synth.allele([*inner[:unit], synth.dupc(), *inner[unit + 1 :]])


def minority(shape: str, af: float, depth: int, seed: int, layout: str) -> list[ReadRecord]:
    """``depth`` reads; ``af`` of the event allele's reads carry dupC at EVENT_UNIT.

    ``two_peak``: SHORT and LONG alleles with ``depth // 2`` reads each, the minority
    inside LONG's peak. ``one_peak``: a homozygous LONG sample. ``af`` 0 is the
    wild-type counterpart (same reads, no minority).
    """
    rng = random.Random(seed)
    alleles = [SHORT, LONG] if layout == "two_peak" else [LONG]
    per = depth // len(alleles)
    out = []
    for a, inner in enumerate(alleles):
        wild, event = synth.allele(inner), carrier(inner)
        n_minor = round(af * per) if inner is LONG else 0
        for i in range(per):
            template = event if i < n_minor else wild
            out.append(read(template, shape, rng, f"m{seed}_{a}_{i}"))
    return out


def artefact(shape: str, share: float, depth: int, seed: int) -> list[ReadRecord]:
    """A homozygous normal whose C unit's C6 run reads C7 in ``share`` of the reads."""
    rng = random.Random(seed)
    wild = synth.allele(ARTEFACT_WT)
    start = base._unit_run(wild, 0)
    edited = wild[:start] + "C" + wild[start:]
    return [
        read(edited if rng.random() < share else wild, shape, rng, f"a{seed}_{i}")
        for i in range(depth)
    ]
