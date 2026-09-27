"""Task 15j: drop within-peak candidate sites explained by a low-accuracy read subset.

On simulated HiFi amplicons, a subset of reads with lower base quality can share a
systematic error at one site, whose minor allele then reaches just above ``het_af_min``
inside one length peak. When each allele already has its own peak, such a site is not
a further allele, yet it left the peak ``unconfirmed_single_site`` (or, with several
such sites linked by the same reads, split it into a third group).

The measure is the read's mean base quality (``SpanRead.mean_q``), which is measured
independently of the consensus. Per-read discordance to the peak draft was rejected:
reads of a second haplotype differ from a majority draft at many sites that are not
candidates, so a true minor haplotype carried by accurate reads would look inaccurate.

A candidate site is dropped only when both hold:

1. its minor carriers have lower mean base quality than its major carriers
   (one-sided Mann-Whitney rank-sum test at ``phase_quality_alpha``), and
2. among the ``phase_quality_keep_frac`` of the site's reads with the highest mean base
   quality, the minor allele fraction is significantly below ``het_af_min``: the
   one-sided upper confidence bound (exact binomial, level ``phase_quality_af_alpha``)
   lies below it. A point estimate would drop a true minor at an allele fraction just
   above ``het_af_min`` by sampling noise alone about half the time whenever condition 1
   holds, for example when insertion stutter of low-quality non-carrier reads enriches
   the minor set of a real homopolymer-run minor in poor reads.

A site whose change is the site-table signature of a dictionary template
(``known_events``: derived by running every template, insertions and deletions, through
the site table; e.g. dupC's C7 -> C8, insG's C7 -> C6, delinsAT's C7 -> C4, dupA's
inserted A) is never dropped (Task 15k): a real dupC minority carried only by
low-quality reads, plus insertion stutter of other low-quality reads at that run, was
otherwise explained away and released to NEGATIVE.

It fails closed: when that high-quality subset has fewer than
``quality_floor_reads(settings)`` reads, when the qualities do not vary, or when any
read of the site table carries no base-quality information (mean Phred 0), the site is
kept.
The caller applies the rule only to a peak of a two-peak model. Reads are never removed
from peak support, allele consensus or event evidence; only the site list changes.
"""

from __future__ import annotations

import math
from typing import Any

from muc_one_span.hybrid.known_events import KnownEventSites, is_known_event_site
from muc_one_span.hybrid.phase_sites import AF_DECIMALS, Site
from muc_one_span.settings import HybridSettings

# Continuity correction of the normal approximation to the rank-sum statistic
# (half a unit of U; part of the test's definition, not a tunable).
CONTINUITY = 0.5
# Reporting precision (significant digits) of a dropped site's p values.
P_DIGITS = 3
# Mean Phred of a read without base-quality information ('!' or missing qualities).
NO_QUALITY = 0.0


def _ranks(values: list[float]) -> tuple[list[float], list[int]]:
    """Mid-ranks (1-based) of ``values`` and the sizes of their tie groups."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    ties = []
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        for k in order[i:j]:
            ranks[k] = (i + j + 1) / 2
        ties.append(j - i)
        i = j
    return ranks, ties


def rank_sum_p_lower(x: list[float], y: list[float]) -> float:
    """One-sided p value that ``x`` tends to be lower than ``y`` (Mann-Whitney U).

    Normal approximation with tie correction and continuity correction. Without
    variance (an empty sample, or every value tied) there is no evidence: 1.0.
    """
    n1, n2 = len(x), len(y)
    n = n1 + n2
    if not n1 or not n2:
        return 1.0
    ranks, ties = _ranks(x + y)
    u = sum(ranks[:n1]) - n1 * (n1 + 1) / 2
    tie_term = sum(t**3 - t for t in ties) / (n * (n - 1))
    var = n1 * n2 / 12 * ((n + 1) - tie_term)
    if var <= 0:
        return 1.0
    z = (u - n1 * n2 / 2 + CONTINUITY) / math.sqrt(var)
    return 0.5 * math.erfc(-z / math.sqrt(2))


def binomial_lower_tail(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p), summed in log space (exact)."""
    if k >= n:
        return 1.0
    log_p, log_q = math.log(p), math.log1p(-p)
    return sum(
        math.exp(
            math.lgamma(n + 1)
            - math.lgamma(i + 1)
            - math.lgamma(n - i + 1)
            + i * log_p
            + (n - i) * log_q
        )
        for i in range(k + 1)
    )


def quality_floor_reads(settings: HybridSettings) -> int:
    """Fewest high-quality reads that can re-test a site (fail closed below it).

    Enough reads to hold ``phase_min_minor_reads`` minor reads at ``het_af_min``, the
    candidate floors themselves.
    """
    return math.ceil(settings.phase_min_minor_reads / settings.het_af_min)


def _explained(
    site: dict[str, Any], feats: list[dict[Site, Any]], quals: list[float], s: HybridSettings
) -> dict[str, Any] | None:
    """The site with its quality evidence when a low-accuracy subset explains it."""
    if any(q <= NO_QUALITY for q in quals):
        return None
    alleles = [(f[site["site"]], q) for f, q in zip(feats, quals, strict=True) if site["site"] in f]
    minor = [q for a, q in alleles if a == site["minor"]]
    major = [q for a, q in alleles if a == site["major"]]
    p = rank_sum_p_lower(minor, major)
    if p >= s.phase_quality_alpha:
        return None
    ranked = sorted(alleles, key=lambda aq: -aq[1])  # stable: input order breaks ties
    kept = ranked[: math.ceil(s.phase_quality_keep_frac * len(ranked))]
    if len(kept) < quality_floor_reads(s):
        return None
    k = sum(a == site["minor"] for a, _q in kept)
    tail = binomial_lower_tail(k, len(kept), s.het_af_min)
    if tail >= s.phase_quality_af_alpha:
        return None
    return {
        **site,
        "af_high_quality": round(k / len(kept), AF_DECIMALS),
        "quality_p": float(f"{p:.{P_DIGITS}g}"),
        "af_bound_p": float(f"{tail:.{P_DIGITS}g}"),
    }


def quality_sites(
    sites: list[dict[str, Any]],
    feats: list[dict[Site, Any]],
    quals: list[float],
    settings: HybridSettings,
    *,
    meta: dict[Site, tuple[str, int]] | None = None,
    insertions: KnownEventSites = KnownEventSites(),
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Split candidate ``sites`` into (kept, dropped as explained by poor reads).

    ``feats`` and ``quals`` are the site table's reads and their mean base qualities, in
    the same order. ``phase_quality_alpha`` 0 keeps every site. A site whose change is
    the site-table signature of a dictionary template (``insertions``, from
    ``known_events.known_event_sites``; run sites need the site table's ``meta`` for
    their base) is always kept (Task 15k).
    """
    kept, dropped = [], []
    for site in sites:
        if is_known_event_site(site, meta or {}, insertions):
            kept.append(site)
            continue
        explained = _explained(site, feats, quals, settings)
        if explained is None:
            kept.append(site)
        else:
            dropped.append(explained)
    return kept, dropped
