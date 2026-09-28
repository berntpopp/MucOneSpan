"""A within-peak homopolymer-run minority above the expected stutter (Task 15l).

A minority haplotype inside a length peak (mosaicism, a third haplotype, a chimera)
that carries a run-length frameshift (dupC: an X unit's C7 run read as C8) at an allele
fraction of 0.15-0.30 stays below every earlier run floor: a candidate site needs
``max(het_af_min, phase_run_bg_multiplier x peer background)`` and the Task 15g safety
tier ``max(het_af_min, phase_run_safety_multiplier x peer background)``, and the 15g
tier applies to single-peak models only. Both compare the raw share of one length with
a fixed multiple of the peers' share, so a minority whose own reads partly stutter back
to the major length never reaches them.

Here each run of an unsplit peak is tested against the stutter it is **expected** to
show. For the run's modal length M (its clean observations) and each other observed
length m, the reads' clean observations (``polish.run_observation``) are explained as a
two-length mixture: per strand, the probability of an observed length under a true M
and under a true m comes from the Task 15f run-length stutter model built from the
peer runs of the same base in this peak (leave-one-out; ``stutter.nearest_profile``
measures, extrapolates or shifts, and ``stutter.event_profile`` applies the
identifiability guard to m). With ``hp_stutter_model = "shift"`` the pre-15f model is
used instead: the major length's own peer profile, moved by m - M for the minor. The weight of m is the stutter-deconvolved minority share;
its one-sided lower confidence bound at ``phase_run_minor_alpha``
(``run_strand.share_lower_bound``) must reach ``phase_run_minor_min_share``. A wild-type
run whose stutter matches its peers has a weight near 0 however often the stuttered
length is seen. Clean observations drop reads whose bounding base changed, so a
neighbouring substitution (``CCCCCC[A]AG`` read as ``CCCCCCCAG``) is not a +1.

A run that passes is only ever a reason to block a negative call (INCONCLUSIVE with the
located site); it never splits a peak and never creates or supports an event.

Every tunable is a validated ``HybridSettings`` field taken from ``settings``.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from muc_one_span.hybrid.evidence import event_allele_fraction
from muc_one_span.hybrid.known_events import KnownEventSites, is_known_event_site
from muc_one_span.hybrid.phase_sites import AF_DECIMALS, Meta, Site
from muc_one_span.hybrid.run_strand import share_lower_bound
from muc_one_span.hybrid.stutter import (
    STRANDS,
    event_profile,
    nearest_profile,
    shift,
    smooth,
)
from muc_one_span.settings import HybridSettings

Clean = list[dict[Site, int]]


def clean_meta(clean: Clean, meta: Meta) -> Meta:
    """Run sites keyed to their modal clean-observed length (draft length if unobserved)."""
    counts: dict[Site, Counter[int]] = {}
    for obs in clean:
        for site, k in obs.items():
            counts.setdefault(site, Counter())[k] += 1
    return {
        site: (base, counts[site].most_common(1)[0][0] if counts.get(site) else length)
        for site, (base, length) in meta.items()
    }


class _Peers:
    """Per-strand clean observations of every run, grouped by (base, modal length)."""

    def __init__(self, clean: Clean, strands: list[str], meta: Meta) -> None:
        self.meta = meta
        self.site: dict[Site, dict[str, Counter[int]]] = {}
        for obs, strand in zip(clean, strands, strict=True):
            for site, k in obs.items():
                self.site.setdefault(site, {}).setdefault(strand, Counter())[k] += 1
        self.cls: dict[tuple[str, int], list[Site]] = {}
        for site, key in meta.items():
            self.cls.setdefault(key, []).append(site)

    def counter(
        self, base: str, length: int, strand: str, exclude: Site
    ) -> tuple[int, Counter[int]]:
        """(peer runs, their observations on ``strand``) of base x length, without ``exclude``."""
        runs = [p for p in self.cls.get((base, length), []) if p != exclude]
        total: Counter[int] = Counter()
        for p in runs:
            total.update(self.site.get(p, {}).get(strand, Counter()))
        return len(runs), total


def _capped(counter: Counter[int], s: HybridSettings) -> Counter[int]:
    out: Counter[int] = Counter()
    for k, n in counter.items():
        out[min(k, s.hp_max_run_len)] += n
    return out


def _profiles(
    peers: _Peers, site: Site, major: int, minor: int, strand: str, s: HybridSettings
) -> tuple[list[float], list[float]]:
    """(minor, major) stutter profiles on one strand from the peers of ``site``."""
    base = peers.meta[site][0]

    def measured(length: int) -> list[float] | None:
        n_runs, c = peers.counter(base, length, strand, site)
        if n_runs < s.hp_stutter_min_class_runs or sum(c.values()) < s.hp_stutter_min_class_reads:
            return None
        return smooth(_capped(c, s), s)

    raw = smooth(_capped(peers.counter(base, major, strand, site)[1], s), s)
    if s.hp_stutter_model == "shift":
        # The pre-15f model (rule 4): the major length's own peers, moved to the minor.
        return shift(raw, minor - major), raw
    major_profile = nearest_profile(measured, major, s) or raw
    return event_profile(measured, major_profile, minor, major, s), major_profile


def _mixture_obs(
    peers: _Peers,
    clean: Clean,
    strands: list[str],
    site: Site,
    alleles: tuple[int, int],
    s: HybridSettings,
) -> list[tuple[float, float]]:
    """(P(obs | minor length), P(obs | major length)) per clean observation of ``site``."""
    major, minor = alleles
    profiles = {st: _profiles(peers, site, major, minor, st, s) for st in {*STRANDS, *strands}}
    obs = []
    for o, strand in zip(clean, strands, strict=True):
        k = o.get(site)
        if k is not None:
            p_minor, p_major = profiles[strand]
            i = min(k, s.hp_max_run_len)
            obs.append((p_minor[i], p_major[i]))
    return obs


def _test(
    peers: _Peers,
    clean: Clean,
    strands: list[str],
    site: Site,
    alleles: tuple[int, int],
    s: HybridSettings,
) -> tuple[float, float]:
    """(stutter-deconvolved minor share, its one-sided lower confidence bound)."""
    obs = _mixture_obs(peers, clean, strands, site, alleles, s)
    share = event_allele_fraction(obs)
    if share < s.phase_run_minor_min_share:
        return share, 0.0  # the bound lies below the share: it cannot pass
    return share, share_lower_bound(obs, s.phase_run_minor_alpha)


def length_aware_share_bound(
    observed: Clean,
    strands: list[str],
    meta: Meta,
    site: Site,
    alleles: tuple[int, int],
    alpha: float,
    s: HybridSettings,
) -> float | None:
    """Lower confidence bound (``alpha``) of the share of a run site's minor length.

    ``observed`` holds each read's run lengths by run site and ``alleles`` is (major,
    minor) length. Each length is convolved with its own Task 15f stutter profile from
    the peer runs of the same base (the tier's model; the single-event gate passes
    the site table's run lengths). None when the minor length is observed in fewer
    than ``phase_min_minor_reads`` reads.
    """
    if sum(o.get(site) == alleles[1] for o in observed) < s.phase_min_minor_reads:
        return None
    peers = _Peers(observed, strands, clean_meta(observed, meta))
    return share_lower_bound(_mixture_obs(peers, observed, strands, site, alleles, s), alpha)


def run_minor_sites(
    clean: Clean,
    strands: list[str],
    meta: Meta,
    settings: HybridSettings,
    known: KnownEventSites,
) -> list[dict[str, Any]]:
    """Run sites whose minority length is significantly above the expected stutter.

    ``clean`` holds the peak's per-read clean run observations (``features(clean=...)``)
    and ``meta`` the site table's runs. For each run and each other observed length
    with at least ``phase_min_minor_reads`` reads within ``phase_run_bg_window`` bases
    of the modal length, the lower confidence bound of the deconvolved share must
    reach ``phase_run_minor_min_share``; with ``phase_run_minor_scope =
    "known_events"`` only lengths that make the change a dictionary template's
    site-table signature (``known_events``) are tested. Largest bound first; each run
    reports its best length.
    """
    if settings.phase_run_minor_scope == "off":
        return []
    modal = clean_meta(clean, meta)
    peers = _Peers(clean, strands, modal)
    out: list[dict[str, Any]] = []
    for site, per in peers.site.items():
        c: Counter[int] = Counter()
        for counter in per.values():
            c.update(counter)
        major = modal[site][1]
        best: tuple[float, float, int] | None = None
        for minor, n in c.items():
            if minor == major or n < settings.phase_min_minor_reads:
                continue
            if abs(minor - major) > settings.phase_run_bg_window:
                continue
            candidate = {"site": site, "major": major, "minor": minor}
            if settings.phase_run_minor_scope == "known_events" and not is_known_event_site(
                candidate, modal, known
            ):
                continue
            share, bound = _test(peers, clean, strands, site, (major, minor), settings)
            if bound >= settings.phase_run_minor_min_share and (best is None or bound > best[1]):
                best = (share, bound, minor)
        if best is None:
            continue
        tot = sum(c.values())
        out.append(
            {
                "site": site,
                "major": major,
                "minor": best[2],
                "af": round(c[best[2]] / tot, AF_DECIMALS),
                "n": tot,
                "share": round(best[0], AF_DECIMALS),
                "share_lower_bound": round(best[1], AF_DECIMALS),
            }
        )
    out.sort(key=lambda site: -site["share_lower_bound"])
    return out
