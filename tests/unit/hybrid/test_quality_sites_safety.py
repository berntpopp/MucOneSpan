"""Task 15j fix round 1: the low-accuracy-subset rule must never release a real minor.

The failure shape (review of the first 15j commit): two length peaks, and inside the
long allele's peak a real minority haplotype carrying dupC (C7 -> C8 in one X unit).
On ONT-like reads, lower-quality reads stutter more; their one-base insertion stutter
turns some non-carrier C7 reads into C8, so the C8 (minor) set is enriched in poor
reads and the rank-sum test can find the minor carriers "lower quality" although the
high-quality carriers read C8 correctly. If the minor allele fraction among the
high-quality reads is then compared with ``het_af_min`` as a point estimate, sampling
noise alone drops a true minor at AF 0.20-0.25 part of the time, and the unresolved
site that blocked a negative call disappears: nothing else blocks NEGATIVE, because
the allele's majority consensus carries no event.

The rule now requires a one-sided upper confidence bound on that high-quality AF to
lie below ``het_af_min``. The shapes here are swept over minor AF 0.20, 0.22 and 0.25
and several seeds, on ONT-like (quality-correlated stutter) and HiFi-like reads; none
may ever be NEGATIVE.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.spans import PHRED_OFFSET, ReadRecord
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_single_event as base

S = DEFAULT_SETTINGS.hybrid
SHORT = ["X"] * 20
LONG = ["X"] * 30
CARRIER = ["X"] * 5 + [synth.dupc()] + ["X"] * 24  # dupC: the X unit's C7 -> C8
N_PER_ALLELE = S.phase_max_site_reads
MINOR_AFS = (0.20, 0.22, 0.25)
SEEDS = tuple(range(6))
# Share of each allele's reads that are low quality.
POOR_FRAC = 0.4
# Share of low-quality non-carrier reads whose one-base insertion stutter at the dupC
# run reads C8 like a carrier. It is site-specific (no other run stutters), which is
# the adversarial case: the peer-run background stays low, so the site is still a
# candidate, while its minor set is strongly enriched in low-quality reads.
POOR_STUTTER_AT_EVENT = 0.3
# Per-read mean Phred ranges by read quality (ONT-like: disjoint; HiFi-like: overlapping).
PROFILES = {
    "ont": {"good": (18, 26), "poor": (8, 17)},
    "hifi": {"good": (30, 40), "poor": (22, 33)},
}


def _read(template: str, rng: random.Random, name: str, q: tuple[int, int]) -> ReadRecord:
    read = synth.reads(template, 1, err=base.ERR, seed=rng.randrange(1 << 30), strand_mix=False)
    seq = read[0].seq if rng.random() < 1 / 2 else rc(read[0].seq)
    return ReadRecord(name, seq, chr(rng.randint(*q) + PHRED_OFFSET) * len(seq))


def _sample(profile: str, af: float, seed: int) -> list[ReadRecord]:
    """Two alleles; ``af`` of the long allele's reads carry dupC (a real minor)."""
    rng = random.Random(seed)
    short, long_, carrier = (synth.allele(a) for a in (SHORT, LONG, CARRIER))
    out = []
    for i in range(N_PER_ALLELE):
        poor = rng.random() < POOR_FRAC
        out.append(_read(short, rng, f"a{i}", PROFILES[profile]["poor" if poor else "good"]))
    n_minor = round(af * N_PER_ALLELE)
    for i in range(N_PER_ALLELE):
        poor = rng.random() < POOR_FRAC
        is_carrier = i < n_minor
        stutter = poor and not is_carrier and rng.random() < POOR_STUTTER_AT_EVENT
        template = carrier if is_carrier or stutter else long_
        out.append(_read(template, rng, f"b{i}", PROFILES[profile]["poor" if poor else "good"]))
    return out


@pytest.mark.parametrize("profile", sorted(PROFILES))
@pytest.mark.parametrize("af", MINOR_AFS)
@pytest.mark.parametrize("seed", SEEDS)
def test_true_within_peak_run_minor_is_never_negative(
    tmp_path: Path, profile: str, af: float, seed: int
) -> None:
    summary, decision = base._run(tmp_path, _sample(profile, af, seed))
    assert decision["state"] != "NO_PATHOGENIC_VARIANT_DETECTED", summary["hybrid"]
