"""Template and read counts that reach a target spanning depth per allele.

Fragment model assumption: MucOneUp ≥ 0.45 fragment sampler starts uniformly
in [-(L-1), S) and is clipped at both ends, where L is the fragment length
and S is the source sequence length. A fragment spans the target VNTR when
its start <= span_start and start + length >= span_end.

PCR slope convention: ``amount.pcr_slope_per_unit`` holds the magnitude of the
(negative) per-unit slope of the log allele ratio for each PCR level; the minor
(shorter) allele share is ``1 / (1 + exp(slope * Δ))``. All tunables come from
`bench_config.AmountConfig`.
"""

from __future__ import annotations

import math
import random

from .bench_config import DEFAULT_BENCH_CONFIG, AmountConfig

MAX_MINOR_SHARE = 0.5  # the minor allele's share by definition never exceeds one half


def pcr_minor_share(
    lengths: tuple[int, int], pcr: str, config: AmountConfig = DEFAULT_BENCH_CONFIG.amount
) -> float:
    """Calculate the expected PCR share of the minor (shorter) allele.

    Args:
        lengths: Tuple of (allele1_length, allele2_length) in bp.
        pcr: PCR level, a key of ``config.pcr_slope_per_unit``.
        config: Amount settings.

    Returns:
        Expected PCR share of the minor (shorter) allele in (0, 0.5].
    """
    delta = abs(lengths[0] - lengths[1])
    return 1.0 / (1.0 + math.exp(config.pcr_slope_per_unit[pcr] * delta))


def capped_minor_share(
    share: float, config: AmountConfig = DEFAULT_BENCH_CONFIG.amount
) -> tuple[float, bool]:
    """Return ``(max(share, config.min_minor_share), capped)`` for sizing template counts.

    Below the floor (strong PCR bias, large length difference) the template count
    explodes and simulation does not finish in bounded time; the minor allele then
    gets below-target depth instead, as in real allelic dropout.
    """
    floor = config.min_minor_share
    return (floor, True) if share < floor else (share, False)


def amplicon_templates(
    target_full_per_allele: int, artefact_rate: float, minor_share: float
) -> int:
    """Calculate template molecules needed to reach target spanning reads per allele.

    The calculation accounts for PCR bias (minor_share) and artefact loss (smear,
    chimera, concatemer), yielding the number of full-length template molecules
    required as input.

    Args:
        target_full_per_allele: Target number of full-length reads for the minor allele.
        artefact_rate: Fraction of molecules lost to artefacts (smear + chimera + concatemer);
                       must be in [0, 1).
        minor_share: Expected PCR share of the minor allele; must be in (0, 0.5].

    Returns:
        Number of template molecules (--coverage) required to reach the target.

    Raises:
        ValueError: If minor_share is not in (0, 0.5] or artefact_rate is not in [0, 1).
    """
    if not 0 < minor_share <= MAX_MINOR_SHARE or not 0 <= artefact_rate < 1:
        raise ValueError(
            f"minor_share in (0, {MAX_MINOR_SHARE}] and artefact_rate in [0, 1) required"
        )
    return math.ceil(target_full_per_allele / (minor_share * (1.0 - artefact_rate)))


def genomic_reads(
    target_spanning: int,
    source_len: int,
    span_start: int,
    span_end: int,
    median: float,
    sigma: float,
    n_hap: int = 2,
    seed: int = 0,
    config: AmountConfig = DEFAULT_BENCH_CONFIG.amount,
) -> int:
    """Calculate genomic reads needed to reach target spanning reads (Monte-Carlo).

    Uses the MucOneUp ≥ 0.45 fragment model: fragment start is uniform in
    [-(L-1), S) and clipped at both ends, where L is fragment length and S
    is the source length. A fragment spans the target VNTR region when
    start <= span_start and start + length >= span_end.

    Args:
        target_spanning: Target number of spanning reads per haplotype.
        source_len: Length of the source sequence (genomic background + VNTR region) in bp.
        span_start: Start position of the target VNTR region (anchor-to-anchor left boundary).
        span_end: End position of the target VNTR region (anchor-to-anchor right boundary).
        median: Median fragment length (lognormal distribution).
        sigma: Shape parameter of the lognormal distribution (scale).
        n_hap: Number of haplotypes (default 2).
        seed: Random seed for reproducibility (default 0).
        config: Amount settings; ``genomic_mc_draws`` Monte-Carlo samples are drawn.

    Returns:
        Number of genomic reads required to reach the target spanning depth per haplotype.

    Raises:
        ValueError: If the fragment length model cannot span the target VNTR region
                    (no spanning reads in the Monte-Carlo sample).
    """
    rng = random.Random(seed)
    draws = config.genomic_mc_draws
    hits = 0
    for _ in range(draws):
        length = max(1, round(rng.lognormvariate(math.log(median), sigma)))
        start = rng.randrange(-(length - 1), source_len)
        hits += start <= span_start and start + length >= span_end
    if hits == 0:
        raise ValueError("fragment length model cannot span the VNTR; lengthen reads or flanks")
    return math.ceil(target_spanning * n_hap * draws / hits)
