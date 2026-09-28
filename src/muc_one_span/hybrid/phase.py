"""S4: split a length peak only when >= min_linked_sites linked, strand-consistent events agree.

Ported from the prototype ``hetsplit.py`` (site table, candidate sites, linkage); the
site table lives in ``phase_sites``. Deviation: groups are formed by a majority vote
over the linked events; the prototype's EM read phasing (``phase_reads``) is not ported.

A site's "minor" allele is chosen per site, so two linked sites can point in opposite
directions. Each linked site is therefore oriented against a reference site (the most
balanced, largest af * (1 - af)) by the sign of phi along the linkage graph before any
read votes. Sites are counted per consensus column, so one unit type that differs from
its neighbour at several non-touching bases yields several linked events (D5).

Every tunable is a validated ``HybridSettings`` field taken from the ``settings``
argument; nothing is defaulted from ``DEFAULT_SETTINGS``.
"""

from __future__ import annotations

import math
import random
from collections import deque
from dataclasses import dataclass, field
from typing import Any

from muc_one_span.hybrid.known_events import KnownEventSites
from muc_one_span.hybrid.phase_quality import guarded_but_explained, quality_sites
from muc_one_span.hybrid.phase_sites import (
    Meta,
    Site,
    candidates,
    events,
    features,
    run_excess_sites,
    strand_biased_sites,
    top_site,
)
from muc_one_span.hybrid.run_minor import Clean, run_minor_sites
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

# Basis of a peak left unsplit with fewer linked events than min_linked_sites.
UNCONFIRMED_SINGLE_SITE = "unconfirmed_single_site"


@dataclass
class PhaseResult:
    """Groups of members and the basis of the (non-)split.

    ``basis`` is "none", "linked_sites", "unconfirmed_single_site" (fewer linked events
    than ``min_linked_sites``) or "unconfirmed_group_size" (a linked split whose smaller
    group is below ``het_min_group``). ``unassigned`` holds members of a split that are
    informative at no linked event or tie between the groups; the caller reassigns them.
    ``run_excess`` holds, for a peak with no candidate site ("none"), the run sites that
    clear only the lower safety floor (``phase_sites.run_excess_sites``), and
    ``strand_biased`` its column and insertion sites refused only for strand bias
    (``phase_sites.strand_biased_sites``); they never split the peak, and the engine
    uses them only to block a negative call. ``quality_associated`` holds the candidate
    sites of an unsplit peak dropped as explained by a low-accuracy read subset
    (``phase_quality.quality_sites``, only when ``quality_filter`` was set).
    ``quality_single_event`` (Task 15k, only when ``single_event_quality`` was set) is
    the unconfirmed result of the same peak without such sites, offered when the sites
    kept form exactly one event; the caller may split on it (``single_event``).
    """

    groups: list[list[SpanRead]]
    basis: str
    sites: list[dict[str, Any]] = field(default_factory=list)
    candidate: dict[str, Any] | None = None
    unassigned: list[SpanRead] = field(default_factory=list)
    run_excess: list[dict[str, Any]] = field(default_factory=list)
    strand_biased: list[dict[str, Any]] = field(default_factory=list)
    run_minor: list[dict[str, Any]] = field(default_factory=list)
    quality_associated: list[dict[str, Any]] = field(default_factory=list)
    quality_single_event: PhaseResult | None = None


def _signed_phi(pairs: list[tuple[int, int]]) -> float:
    """Signed phi correlation of two binary indicators over the reads informative for both."""
    a = sum(1 for x, y in pairs if x and y)
    b = sum(1 for x, y in pairs if x and not y)
    c = sum(1 for x, y in pairs if not x and y)
    d = len(pairs) - a - b - c
    den = math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
    return (a * d - b * c) / den if den else 0.0


def _indicator(f: dict[Site, Any], s: dict[str, Any]) -> int | None:
    """1 for the site's minor allele, 0 for its major, None when uninformative."""
    a = f.get(s["site"])
    return 1 if a == s["minor"] else 0 if a == s["major"] else None


def _pairs(vec: list[list[int | None]], i: int, j: int) -> list[tuple[int, int]]:
    out = []
    for row in vec:
        x, y = row[i], row[j]
        if x is not None and y is not None:
            out.append((x, y))
    return out


def _components(adj: dict[int, dict[int, int]]) -> list[dict[int, int]]:
    """Connected components, each as {site index: orientation relative to its first site}."""
    seen: set[int] = set()
    comps = []
    for start in adj:
        if start in seen:
            continue
        ori, queue = {start: 1}, deque([start])
        seen.add(start)
        while queue:
            k = queue.popleft()
            for m, sign in adj[k].items():
                if m not in seen:
                    seen.add(m)
                    ori[m] = ori[k] * sign
                    queue.append(m)
        comps.append(ori)
    return comps


def _linked(
    feats: list[dict[Site, Any]],
    sites: list[dict[str, Any]],
    meta: Meta,
    settings: HybridSettings,
) -> dict[int, int]:
    """Component of phi-linked sites with the most events, oriented to its top site.

    Returns {site index: +1 or -1}; -1 means the site's minor allele travels with the
    reference site's major allele, so its indicator is flipped before voting.
    """
    vec = [[_indicator(f, s) for s in sites] for f in feats]
    adj: dict[int, dict[int, int]] = {i: {} for i in range(len(sites))}
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            pairs = _pairs(vec, i, j)
            if len(pairs) < settings.phase_min_pair_reads:
                continue
            phi = _signed_phi(pairs)
            if abs(phi) >= settings.link_phi_min:
                adj[i][j] = adj[j][i] = -1 if phi < 0 else 1
    best: dict[int, int] = {}
    n_best = 0
    for comp in _components(adj):
        n_events = len(events([sites[k] for k in comp], meta))
        if n_events > n_best:
            best, n_best = comp, n_events
    ref = max(best, key=lambda k: sites[k]["af"] * (1 - sites[k]["af"]))
    return {k: o * best[ref] for k, o in best.items()}


def _majority(votes: list[int]) -> int | None:
    """1 or 0 by strict majority; None when empty or tied."""
    ones = sum(votes)
    zeros = len(votes) - ones
    return None if ones == zeros else int(ones > zeros)


def _vote(
    f: dict[Site, Any], evs: list[list[int]], sites: list[dict[str, Any]], ori: list[int]
) -> int | None:
    """Group of one read: majority over events of each event's oriented site majority."""
    event_votes = []
    for ev in evs:
        known = []
        for k in ev:
            ind = _indicator(f, sites[k])
            if ind is not None:
                known.append(ind if ori[k] > 0 else 1 - ind)
        vote = _majority(known)
        if vote is not None:
            event_votes.append(vote)
    return _majority(event_votes)


def split_by_linked_sites(
    cons: str,
    members: list[SpanRead],
    settings: HybridSettings,
    rng: random.Random,
    *,
    quality_filter: bool = False,
    single_event_quality: bool = False,
    insertions: KnownEventSites | None = None,
) -> PhaseResult:
    """Return one group (no split) unless >= min_linked_sites linked events support two.

    Sites are found on at most ``phase_max_site_reads`` members (sampled with ``rng``);
    every member is then placed by its oriented event votes. Members with no informative
    event or a tied vote go to ``unassigned``.

    With ``quality_filter`` (the caller sets it for a peak of a two-peak model only), a
    peak left unsplit is re-examined without the candidate sites explained by a
    low-accuracy read subset (``phase_quality.quality_sites``; a site adding one of the
    dictionary event signatures ``insertions`` is never dropped), and that result is used
    when it also leaves the peak as one group. A linked split is never undone and no
    split is created, so the members of every group, and hence every allele consensus,
    are the same as without the filter; only the unsplit peak's basis can change.

    With ``single_event_quality`` (Task 15k; the caller sets it for a single-peak
    model only), a peak left ``unconfirmed_single_site`` by more than one event is
    re-examined without the sites the same test explains. When the sites kept form
    exactly one event, that unconfirmed result is attached as ``quality_single_event``
    (with the dropped sites in its ``quality_associated``); the returned result itself
    is unchanged, so the caller decides whether a single-event split is made on it.
    """
    cap = settings.phase_max_site_reads
    sample = members if len(members) <= cap else rng.sample(members, cap)
    clean: Clean = []
    feats, meta = features(cons, [m.seq for m in sample], settings, clean=clean)
    strands = [m.strand for m in sample]
    sites = candidates(feats, strands, meta, settings)
    known = insertions or KnownEventSites()
    table = _SiteTable(cons, members, sample, feats, meta, strands, known, clean)
    result = _split(table, sites, settings)
    if single_event_quality and result.basis == UNCONFIRMED_SINGLE_SITE:
        result.quality_single_event = _single_event_alternative(table, sites, settings)
    if not quality_filter or len(result.groups) > 1:
        return result
    kept, dropped = table.quality(sites, settings)
    if not dropped:
        return result
    filtered = _split(table, kept, settings)
    if len(filtered.groups) > 1:
        return result
    filtered.quality_associated = dropped
    return filtered


@dataclass
class _SiteTable:
    """The phase site table of one peak: its members, the sampled reads and their sites."""

    cons: str
    members: list[SpanRead]
    sample: list[SpanRead]
    feats: list[dict[Site, Any]]
    meta: Meta
    strands: list[str]
    insertions: KnownEventSites
    clean: Clean

    def tier_reads(self, settings: HybridSettings) -> tuple[Clean, list[str]]:
        """Clean run observations and strands for the run-minority tier (Task 15l).

        At most ``phase_run_minor_max_reads`` members: the site table's own sample when
        it already holds that many (or every member), else a fresh sample drawn with
        ``random.Random(settings.seed)``, independent of the site table.
        """
        cap = settings.phase_run_minor_max_reads
        if len(self.sample) >= min(cap, len(self.members)):
            return self.clean[:cap], self.strands[:cap]
        reads = random.Random(settings.seed).sample(self.members, min(cap, len(self.members)))
        clean: Clean = []
        features(self.cons, [m.seq for m in reads], settings, clean=clean)
        return clean, [m.strand for m in reads]

    def quality(
        self, sites: list[dict[str, Any]], settings: HybridSettings
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        """``phase_quality.quality_sites`` on this table's sampled reads."""
        quals = [m.mean_q for m in self.sample]
        return quality_sites(
            sites, self.feats, quals, settings, meta=self.meta, insertions=self.insertions
        )

    def guard_only(
        self, sites: list[dict[str, Any]], settings: HybridSettings
    ) -> list[dict[str, Any]]:
        """``phase_quality.guarded_but_explained`` on this table's sampled reads."""
        quals = [m.mean_q for m in self.sample]
        return guarded_but_explained(
            sites, self.feats, quals, settings, meta=self.meta, insertions=self.insertions
        )


def _single_event_alternative(
    t: _SiteTable, sites: list[dict[str, Any]], settings: HybridSettings
) -> PhaseResult | None:
    """The unconfirmed one-event result without low-accuracy sites, or None.

    Only for a peak whose candidate sites form more than one event: the sites the
    15j test explains (``phase_quality.quality_sites``) are dropped, and the result
    is offered only when the sites kept form exactly one event, and never when a kept
    site is protected only by the known-event guard although poor reads explain it
    (fail closed: that event may be an artefact sharing a template's signature).
    """
    if len(events(sites, t.meta)) == 1:
        return None  # one event already: the ordinary single-event split applies
    kept, dropped = t.quality(sites, settings)
    if not dropped or len(events(kept, t.meta)) != 1:
        return None
    if t.guard_only(kept, settings):
        return None
    alternative = _split(t, kept, settings)
    if alternative.basis != UNCONFIRMED_SINGLE_SITE:
        return None
    alternative.quality_associated = dropped
    return alternative


def _split(t: _SiteTable, sites: list[dict[str, Any]], settings: HybridSettings) -> PhaseResult:
    """The phase result of one peak from its candidate ``sites`` (see split_by_linked_sites)."""
    members, feats, meta = t.members, t.feats, t.meta
    if not sites:
        return PhaseResult(
            [members],
            "none",
            run_excess=run_excess_sites(feats, meta, settings),
            strand_biased=strand_biased_sites(feats, t.strands, meta, settings),
            run_minor=run_minor_sites(*t.tier_reads(settings), meta, settings, t.insertions),
        )
    ori = _linked(feats, sites, meta, settings) if len(sites) > 1 else {0: 1}
    order = sorted(ori)
    linked = [sites[k] for k in order]
    evs = events(linked, meta)
    if len(evs) < settings.min_linked_sites:
        return PhaseResult([members], UNCONFIRMED_SINGLE_SITE, sites, candidate=top_site(sites))
    orient = [ori[k] for k in order]
    all_feats = (
        feats if t.sample is members else features(t.cons, [m.seq for m in members], settings)[0]
    )
    groups: list[list[SpanRead]] = [[], []]
    unassigned: list[SpanRead] = []
    for m, f in zip(members, all_feats, strict=True):
        vote = _vote(f, evs, linked, orient)
        (unassigned if vote is None else groups[vote]).append(m)
    if min(len(g) for g in groups) < settings.het_min_group * len(members):
        return PhaseResult([members], "unconfirmed_group_size", linked, candidate=top_site(linked))
    return PhaseResult(groups, "linked_sites", linked, unassigned=unassigned)
