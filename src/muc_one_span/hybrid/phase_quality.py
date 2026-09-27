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

A site whose minor allele is a dictionary insertion in its context (``known_insertions``:
a run of a template's base and parent-unit length gaining the template's copies, e.g.
dupC's C7 -> C8, or an insertion slot adding a template's sequence) is never dropped
(Task 15k): a real dupC minority carried only by low-quality reads, plus insertion
stutter of other low-quality reads at that run, was otherwise explained away and
released to NEGATIVE.

It fails closed: when that high-quality subset has fewer than
``quality_floor_reads(settings)`` reads, when the qualities do not vary, or when any
read of the site table carries no base-quality information (mean Phred 0), the site is
kept.
The caller applies the rule only to a peak of a two-peak model. Reads are never removed
from peak support, allele consensus or event evidence; only the site list changes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from muc_one_span.config import RepeatDictionary
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


@dataclass(frozen=True)
class KnownInsertions:
    """Insertion events of the repeat dictionary, as the site table sees them.

    ``runs`` holds (base, run length, inserted copies) for a homopolymer insertion
    that lengthens a run of its own base in a parent unit (dupC: ("C", 7, 1));
    ``inserted`` holds every inserted sequence, for insertion-slot sites.
    """

    runs: frozenset[tuple[str, int, int]] = frozenset()
    inserted: frozenset[str] = frozenset()


def _run_at(parent: str, point: int, base: str) -> int:
    """Length of the run of ``base`` in ``parent`` touching the gap before ``point``."""
    lo = point
    while lo > 0 and parent[lo - 1] == base:
        lo -= 1
    hi = point
    while hi < len(parent) and parent[hi] == base:
        hi += 1
    return hi - lo


def known_insertions(rd: RepeatDictionary) -> KnownInsertions:
    """The dictionary's insertion templates, with each parent unit's run context."""
    runs: set[tuple[str, int, int]] = set()
    inserted: set[str] = set()
    for template in rd.mutations.values():
        for change in template.get("changes") or []:
            seq = str(change.get("sequence") or "").upper()
            if change.get("type") != "insert" or not seq:
                continue
            inserted.add(seq)
            if len(set(seq)) != 1:
                continue
            point = int(change["start"]) - 1  # dictionary coordinates are 1-based
            for unit in template.get("allowed_repeats") or []:
                parent = rd.repeats.get(unit, "").upper()
                length = _run_at(parent, point, seq[0]) if 0 <= point <= len(parent) else 0
                if length:
                    runs.add((seq[0], length, len(seq)))
    return KnownInsertions(frozenset(runs), frozenset(inserted))


def inserts_known_event(
    site: dict[str, Any], meta: dict[Site, tuple[str, int]], known: KnownInsertions
) -> bool:
    """True when the site's minor allele is a dictionary insertion in its context.

    A run site matches when its base, major length and gained copies are those of a
    template's run (dupC: a C7 run read as C8); an insertion slot matches when its
    minor adds exactly a template's inserted sequence. Column sites never match.
    """
    kind = site["site"][0]
    if kind == "run" and site["site"] in meta and site["minor"] > site["major"]:
        base = meta[site["site"]][0]
        return (base, site["major"], site["minor"] - site["major"]) in known.runs
    if kind == "ins" and str(site["minor"]).startswith(str(site["major"])):
        return str(site["minor"])[len(str(site["major"])) :] in known.inserted
    return False


def quality_sites(
    sites: list[dict[str, Any]],
    feats: list[dict[Site, Any]],
    quals: list[float],
    settings: HybridSettings,
    *,
    meta: dict[Site, tuple[str, int]] | None = None,
    insertions: KnownInsertions = KnownInsertions(),
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Split candidate ``sites`` into (kept, dropped as explained by poor reads).

    ``feats`` and ``quals`` are the site table's reads and their mean base qualities, in
    the same order. ``phase_quality_alpha`` 0 keeps every site. A site whose minor
    allele is a dictionary insertion in its context (``insertions``, from
    ``known_insertions``; run sites need the site table's ``meta`` for their base) is
    always kept (Task 15k).
    """
    kept, dropped = [], []
    for site in sites:
        if inserts_known_event(site, meta or {}, insertions):
            kept.append(site)
            continue
        explained = _explained(site, feats, quals, settings)
        if explained is None:
            kept.append(site)
        else:
            dropped.append(explained)
    return kept, dropped
