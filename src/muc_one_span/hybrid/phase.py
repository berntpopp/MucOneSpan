"""S4: split a length peak only when >= min_linked_sites read-linked difference sites agree.

Ported from the prototype ``hetsplit.py`` (site table, candidate sites, linkage).
Deviation: groups are formed by a majority vote over the linked sites; the prototype's
EM read phasing (``phase_reads``) is not ported. Sites are counted per consensus
column, so one unit type that differs from its neighbour at several bases yields
several perfectly linked sites.

Every tunable (site-read cap, homopolymer-run site length and background window, minor
alleles per site and their read floor, run background multiplier, gap AF factor and the
pairwise-linkage read floor) is a validated ``HybridSettings`` field taken from the
``settings`` argument; nothing is defaulted from ``DEFAULT_SETTINGS``.
"""

from __future__ import annotations

import math
import random
import statistics
from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from muc_one_span.hybrid.align import global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

AF_DECIMALS = 3  # reporting precision of a site's allele fraction (format, not a tunable)
GAP = "-"  # deletion symbol in projected columns (``align.Columns``)
Site = tuple[str, int]


@dataclass
class PhaseResult:
    """Groups of members and the basis of the (non-)split."""

    groups: list[list[SpanRead]]
    basis: str  # "none", "linked_sites" or "unconfirmed_single_site"
    sites: list[dict[str, Any]] = field(default_factory=list)
    candidate: dict[str, Any] | None = None


def _features(
    cons: str, reads: list[str], settings: HybridSettings
) -> tuple[list[dict[Site, Any]], dict[Site, tuple[str, int]]]:
    """Per-read alleles at every column, insertion slot and homopolymer-run site."""
    runs = _runs(cons, settings.phase_run_min_len)
    in_run = {i for s, e, _ in runs for i in range(s, e)}
    feats = []
    for read in reads:
        proj = global_columns(read, cons)
        f: dict[Site, Any] = {}
        for pos, base in enumerate(proj.cols):
            if pos not in in_run:
                f[("col", pos)] = base
                if pos - 1 not in in_run:
                    f[("ins", pos)] = proj.ins.get(pos, "")
        for s, e, b in runs:
            f[("run", s)] = read_run_length(read, proj.t2q, s, e, b)
        feats.append(f)
    return feats, {("run", s): (b, e - s) for s, e, b in runs}


def _run_background(
    counts: dict[Site, Counter[Any]], meta: dict[Site, tuple[str, int]], window: int
) -> dict[tuple[str, int, int], float]:
    """Median fraction, over runs of the same base and length, of reads showing length x."""
    fracs: dict[tuple[str, int, int], list[float]] = {}
    for site, (base, length) in meta.items():
        c = counts.get(site, Counter())
        tot = sum(c.values()) or 1
        for x in set(range(max(0, length - window), length + window + 1)) - {length}:
            fracs.setdefault((base, length, x), []).append(c.get(x, 0) / tot)
    return {k: statistics.median(v) for k, v in fracs.items()}


def _candidates(
    feats: list[dict[Site, Any]], meta: dict[Site, tuple[str, int]], settings: HybridSettings
) -> list[dict[str, Any]]:
    """Sites whose best minor allele clears het_af_min (and the run/gap noise rules)."""
    counts: dict[Site, Counter[Any]] = {}
    for f in feats:
        for site, allele in f.items():
            counts.setdefault(site, Counter())[allele] += 1
    bg = _run_background(counts, meta, settings.phase_run_bg_window)
    mult, af_min = settings.phase_run_bg_multiplier, settings.het_af_min
    out = []
    for site, c in counts.items():
        tot = sum(c.values())
        major, _ = c.most_common(1)[0]
        minors = [
            (a, n)
            for a, n in c.most_common(settings.phase_max_minor_alleles + 1)
            if a != major and n >= settings.phase_min_minor_reads
        ]
        if not minors:
            continue
        if site[0] == "run":
            base, length = meta[site]
            _score, minor, n_minor = max(
                (n / tot - mult * bg.get((base, length, a), 0.0), a, n) for a, n in minors
            )
            if n_minor / tot < max(af_min, mult * bg.get((base, length, minor), 0.0)):
                continue
        else:
            minor, n_minor = minors[0]
            af = n_minor / tot
            if af < af_min or (
                GAP in (major, minor) and af < settings.phase_gap_af_factor * af_min
            ):
                continue
        out.append(
            {
                "site": site,
                "major": major,
                "minor": minor,
                "af": round(n_minor / tot, AF_DECIMALS),
                "n": tot,
            }
        )
    return out


def _phi(pairs: list[tuple[int, int]]) -> float:
    """|phi| correlation of two binary indicators over the reads informative for both."""
    a = sum(1 for x, y in pairs if x and y)
    b = sum(1 for x, y in pairs if x and not y)
    c = sum(1 for x, y in pairs if not x and y)
    d = len(pairs) - a - b - c
    den = math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
    return abs(a * d - b * c) / den if den else 0.0


def _indicator(f: dict[Site, Any], s: dict[str, Any]) -> int | None:
    """1 for the site's minor allele, 0 for its major, None when uninformative."""
    a = f.get(s["site"])
    return 1 if a == s["minor"] else 0 if a == s["major"] else None


def _linked(
    feats: list[dict[Site, Any]], sites: list[dict[str, Any]], settings: HybridSettings
) -> list[dict[str, Any]]:
    """Largest connected component of sites whose read indicators correlate (|phi|)."""
    vec = [[_indicator(f, s) for s in sites] for f in feats]
    n = len(sites)
    adj: dict[int, set[int]] = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            pairs = [(x[i], x[j]) for x in vec if x[i] is not None and x[j] is not None]
            if len(pairs) >= settings.phase_min_pair_reads and (
                _phi(pairs) >= settings.link_phi_min  # type: ignore[arg-type]
            ):
                adj[i].add(j)
                adj[j].add(i)
    seen: set[int] = set()
    best: list[int] = []
    for i in range(n):
        if i in seen:
            continue
        comp, stack = [], [i]
        seen.add(i)
        while stack:
            k = stack.pop()
            comp.append(k)
            for m in adj[k] - seen:
                seen.add(m)
                stack.append(m)
        best = comp if len(comp) > len(best) else best
    return [sites[i] for i in sorted(best)]


def split_by_linked_sites(
    cons: str, members: list[SpanRead], settings: HybridSettings, rng: random.Random
) -> PhaseResult:
    """Return one group (no split) unless >= min_linked_sites linked sites support two.

    Sites are found on at most ``phase_max_site_reads`` members (sampled with ``rng``);
    every member is then placed by a strict majority of its informative linked sites.
    Members informative at no linked site are left out of both groups. A split whose
    smaller group is below ``het_min_group`` of the members is reported unconfirmed.
    """
    cap = settings.phase_max_site_reads
    sample = members if len(members) <= cap else rng.sample(members, cap)
    feats, meta = _features(cons, [m.seq for m in sample], settings)
    sites = _candidates(feats, meta, settings)
    if not sites:
        return PhaseResult([members], "none")
    linked = _linked(feats, sites, settings) if len(sites) > 1 else sites[:1]
    if len(linked) < settings.min_linked_sites:
        top = max(sites, key=lambda s: s["af"] * (1 - s["af"]))
        return PhaseResult([members], "unconfirmed_single_site", sites, candidate=top)
    all_feats = (
        feats if sample is members else _features(cons, [m.seq for m in members], settings)[0]
    )
    groups: list[list[SpanRead]] = [[], []]
    for m, f in zip(members, all_feats, strict=True):
        known = [v for v in (_indicator(f, s) for s in linked) if v is not None]
        if known:
            groups[int(sum(known) * 2 > len(known))].append(m)
    if min(len(g) for g in groups) < settings.het_min_group * len(members):
        return PhaseResult([members], "unconfirmed_single_site", linked, candidate=linked[0])
    return PhaseResult(groups, "linked_sites", linked)
