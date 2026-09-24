import math

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, AmountConfig
from muc_one_span.benchsim.depth import (
    amplicon_templates,
    capped_minor_share,
    genomic_reads,
    pcr_minor_share,
)

AMOUNT = DEFAULT_BENCH_CONFIG.amount


def test_minor_share() -> None:
    assert pcr_minor_share((40, 40), "calibrated") == 0.5
    assert pcr_minor_share((40, 60), "none") == 0.5
    slope = AMOUNT.pcr_slope_per_unit["calibrated"]
    assert pcr_minor_share((40, 60), "calibrated") == pytest.approx(1 / (1 + math.exp(slope * 20)))
    flat = AmountConfig(pcr_slope_per_unit={"calibrated": 0.0, "strong": 0.0, "none": 0.0})
    assert pcr_minor_share((40, 60), "strong", flat) == 0.5
    assert pcr_minor_share((40, 60), "strong") < pcr_minor_share((40, 60), "calibrated")


def test_amplicon_templates_reach_target() -> None:
    n = amplicon_templates(30, 0.287, 0.25)
    assert n * 0.25 * (1 - 0.287) >= 30 > (n - 1) * 0.25 * (1 - 0.287)


def test_amplicon_templates_rejects_invalid_minor_share() -> None:
    with pytest.raises(ValueError, match="minor_share in"):
        amplicon_templates(30, 0.287, 0.0)
    with pytest.raises(ValueError, match="minor_share in"):
        amplicon_templates(30, 0.287, 0.6)


def test_genomic_reads_hits_expected_spanning() -> None:
    import random

    n = genomic_reads(20, 26000, 10000, 16000, 6000, 0.5, n_hap=2, seed=1)
    rng, spans = random.Random(2), 0
    for _ in range(n):
        length = max(1, round(rng.lognormvariate(math.log(6000), 0.5)))
        start = rng.randrange(-(length - 1), 26000)
        spans += start <= 10000 and start + length >= 16000
    assert spans / 2 == pytest.approx(20, rel=0.3)


def test_genomic_reads_rejects_impossible_span() -> None:
    with pytest.raises(ValueError, match="cannot span"):
        genomic_reads(10, 30000, 0, 30000, 500, 0.1, seed=1)


def test_capped_minor_share_floor() -> None:
    floor = AMOUNT.min_minor_share
    assert capped_minor_share(floor / 100) == (floor, True)
    assert capped_minor_share(floor * 2) == (floor * 2, False)
    lower = AmountConfig(min_minor_share=floor / 10)
    assert capped_minor_share(floor / 100, lower) == (floor / 10, True)


def test_genomic_draws_are_configured() -> None:
    few = AmountConfig(genomic_mc_draws=10)
    n = genomic_reads(20, 26000, 10000, 16000, 6000, 0.5, seed=1, config=few)
    assert n > 0
