import math

import pytest

from muc_one_span.benchsim.depth import amplicon_templates, genomic_reads, pcr_minor_share


def test_minor_share() -> None:
    assert pcr_minor_share((40, 40), "calibrated") == 0.5
    assert pcr_minor_share((40, 60), "none") == 0.5
    assert pcr_minor_share((40, 60), "calibrated") == pytest.approx(1 / (1 + math.exp(0.056 * 20)))
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
    from muc_one_span.benchsim.depth import MIN_MINOR_SHARE, capped_minor_share

    assert MIN_MINOR_SHARE == 0.05
    assert capped_minor_share(1e-4) == (0.05, True)
    assert capped_minor_share(0.3) == (0.3, False)
