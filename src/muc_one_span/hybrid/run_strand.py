"""Stutter-aware strand-bias test and read classification at homopolymer-run sites.

A run-length site cannot use the exact Fisher allele x strand test of a column site:
nanopore homopolymer stutter is strand-asymmetric (a C run is read with many more
length errors on one strand than on the other), so a true heterozygous run looks
"depleted" on the noisier strand although both strands carry the two alleles in the
same proportion. Here each observed run length is explained, per strand, by the
length-error profile measured at the *other* runs of the same base in the site table
(leave-one-out; same consensus length first, any length of that base otherwise, as in
``phase_sites._run_background``). The allele weight ``f`` of the minor length in the
two-length mixture is estimated per strand and jointly (``evidence.event_allele_fraction``);
the joint weight must reach ``het_af_min`` (the caller's floor), and a likelihood-ratio test (one degree of freedom) of "one weight for both strands"
against "a weight per strand" gives the strand-bias p value. A strand-specific
systematic error at one site (present on one strand only) still fails; strand-specific
stutter that every run of that base shares is absorbed by the profiles.

Every tunable is a validated ``HybridSettings`` field (``hp_max_run_len``,
``hp_background_pseudocount``); nothing is defaulted from ``DEFAULT_SETTINGS``.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any

from muc_one_span.hybrid.evidence import event_allele_fraction
from muc_one_span.settings import HybridSettings

Site = tuple[str, int]
Meta = dict[Site, tuple[str, int]]
# strand -> smoothed probability of each length error (observed - true), index + cap.
ShiftProfile = dict[str, list[float]]


def shift_profiles(
    feats: list[dict[Site, Any]],
    strands: list[str],
    meta: Meta,
    site: Site,
    length: int,
    settings: HybridSettings,
) -> ShiftProfile:
    """Per-strand length-error profile of runs of the site's base with true ``length``.

    Errors beyond ``hp_max_run_len`` are capped; each profile adds
    ``hp_background_pseudocount`` per error value, so a strand without peer
    observations is uniform (uninformative) rather than missing.
    """
    base = meta[site][0]
    others = [p for p, (b, _n) in meta.items() if p != site and b == base]
    peers = [p for p in others if meta[p][1] == length] or others
    cap = settings.hp_max_run_len
    per: dict[str, Counter[int]] = {st: Counter() for st in strands}
    for f, strand in zip(feats, strands, strict=True):
        for p in peers:
            observed = f.get(p)
            if observed is not None:
                per[strand][max(-cap, min(cap, observed - meta[p][1]))] += 1
    size = 2 * cap + 1
    pseudo = settings.hp_background_pseudocount
    out = {}
    for strand, c in per.items():
        tot = sum(c.values())
        out[strand] = [(c.get(i - cap, 0) + pseudo) / (tot + pseudo * size) for i in range(size)]
    return out


def length_prob(profile: ShiftProfile, strand: str, observed: int, true_len: int) -> float:
    """P(observed run length | true length) on ``strand`` from its error profile."""
    p = profile[strand]
    cap = len(p) // 2
    return p[max(-cap, min(cap, observed - true_len)) + cap]


def _loglik(obs: list[tuple[float, float]], f: float) -> float:
    return sum(math.log(f * p1 + (1 - f) * p0) for p1, p0 in obs)


def run_mixture(
    feats: list[dict[Site, Any]],
    strands: list[str],
    meta: Meta,
    site: Site,
    alleles: tuple[int, int],
    settings: HybridSettings,
) -> tuple[float, float]:
    """(minor-length weight, strand-bias p value) of a run site's two-length mixture.

    ``alleles`` is (major length, minor length). The weight is the maximum-likelihood
    share of the minor length over all reads once each strand's stutter is modelled, so
    a draft run of the wrong length, or stutter that the other runs share, gives a
    weight near 0 however often the stuttered length is observed. The p value tests
    "one weight for both strands" against "a weight per strand" by the likelihood-ratio
    statistic ``2 * (sum_s max log L_s(f_s) - max log L(f))`` against a chi-square
    distribution with one degree of freedom; with reads on one strand only there is
    nothing to compare and the p value is 1.
    """
    major, minor = alleles
    p_minor = shift_profiles(feats, strands, meta, site, minor, settings)
    p_major = shift_profiles(feats, strands, meta, site, major, settings)
    obs: dict[str, list[tuple[float, float]]] = {}
    for f, strand in zip(feats, strands, strict=True):
        k = f.get(site)
        if k is not None:
            obs.setdefault(strand, []).append(
                (length_prob(p_minor, strand, k, minor), length_prob(p_major, strand, k, major))
            )
    pooled = [o for v in obs.values() for o in v]
    weight = event_allele_fraction(pooled)
    if len(obs) < 2:
        return weight, 1.0
    split = sum(_loglik(v, event_allele_fraction(v)) for v in obs.values())
    stat = max(0.0, 2 * (split - _loglik(pooled, weight)))
    # Chi-square survival function with one degree of freedom (the extra strand weight).
    return weight, math.erfc(math.sqrt(stat / 2))
