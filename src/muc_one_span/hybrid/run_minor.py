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
and under a true m comes from the Task 15f run-length stutter model built from the peer runs
of the same base in this peak (leave-one-out; ``stutter.nearest_profile`` measures,
extrapolates or shifts, and ``stutter.event_profile`` applies the identifiability guard to
m). With ``hp_stutter_model = "shift"`` the pre-15f model is used instead: the major
length's own peer profile, moved by m - M for the minor. The weight of m is the
stutter-deconvolved minority share; its one-sided lower confidence bound at
``phase_run_minor_alpha`` (``run_strand.share_lower_bound``) must reach
``phase_run_minor_min_share``. A wild-type run whose stutter matches its peers has a weight
near 0 however often the stuttered length is seen. Clean observations drop reads whose
bounding base changed, so a neighbouring substitution (``CCCCCC[A]AG`` read as
``CCCCCCCAG``) is not a +1.

**Within-run events.** insG, insG_pos58 and delinsAT put another base inside an X
unit's C7 run, so a carrier read never observes that run cleanly and the clean test
cannot see such a minority. The site table reads a carrier as the run shortened to its
longest stretch (C7 -> C6, C5 or C4), the template's run signature (``known_events``).
With ``phase_run_minor_in_run``, each run's bounded observations (both bounding bases
kept; ``features(impure=...)`` holds the impure ones) are also tested for every
signature length m of the run's (base, modal length): per strand, a read shows m inside
the run with probability b at the wild type (the rate at the peer runs of the same base
and length, leave-one-out, with ``hp_background_pseudocount``) and c at a carrier (the
Task 15f stutter probability that a run of length m reads m). The weight of the carrier
component and its one-sided bound are computed as above against the same floor.

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
# Observation kinds of a tier site (report labels): a clean run length, or a bounded
# but impure run whose longest stretch is a known-event signature.
CLEAN, IN_RUN = "clean", "in_run"
# An observation shows the signature length inside the run or not (structural).
SIGNATURE_OUTCOMES = 2


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
    return _bound(_mixture_obs(peers, clean, strands, site, alleles, s), s)


class _InRun:
    """Per-strand bounded observations and impure lengths of every run (``features``)."""

    def __init__(self, clean: Clean, impure: Clean, strands: list[str]) -> None:
        self.bounded: dict[Site, Counter[str]] = {}
        self.lengths: dict[Site, dict[str, Counter[int]]] = {}
        for c, i, strand in zip(clean, impure, strands, strict=True):
            for site in (*c, *i):
                self.bounded.setdefault(site, Counter())[strand] += 1
            for site, k in i.items():
                self.lengths.setdefault(site, {}).setdefault(strand, Counter())[k] += 1

    def rate(self, runs: list[Site], strand: str, length: int, s: HybridSettings) -> float:
        """Smoothed share of ``runs``' bounded observations impure at ``length``."""
        n = sum(self.bounded.get(p, Counter())[strand] for p in runs)
        k = sum(self.lengths.get(p, {}).get(strand, Counter())[length] for p in runs)
        pseudo = s.hp_background_pseudocount
        return (k + pseudo) / (n + SIGNATURE_OUTCOMES * pseudo)


def _in_run_obs(
    peers: _Peers,
    tally: _InRun,
    observed: tuple[Clean, Clean, list[str]],
    site: Site,
    alleles: tuple[int, int],
    s: HybridSettings,
) -> list[tuple[float, float]]:
    """(P(obs | carrier), P(obs | wild type)) per bounded observation of ``site``."""
    clean, impure, strands = observed
    major, minor = alleles
    base = peers.meta[site][0]
    others = [p for p, (b, _n) in peers.meta.items() if p != site and b == base]
    runs = [p for p in others if peers.meta[p][1] == major] or others
    rates = {}
    for strand in set(strands):
        carrier = _profiles(peers, site, major, minor, strand, s)[0]
        rates[strand] = (
            carrier[min(minor, s.hp_max_run_len)],
            tally.rate(runs, strand, minor, s),
        )
    obs = []
    for c, i, strand in zip(clean, impure, strands, strict=True):
        if site in c or site in i:
            p_carrier, p_wild = rates[strand]
            shows = i.get(site) == minor
            obs.append((p_carrier, p_wild) if shows else (1 - p_carrier, 1 - p_wild))
    return obs


def _bound(obs: list[tuple[float, float]], s: HybridSettings) -> tuple[float, float]:
    """(mixture weight, its one-sided lower bound); 0 bound when the weight is below the floor."""
    share = event_allele_fraction(obs)
    if share < s.phase_run_minor_min_share:
        return share, 0.0  # the bound lies below the share: it cannot pass
    return share, share_lower_bound(obs, s.phase_run_minor_alpha)


def _in_run_sites(
    peers: _Peers,
    observed: tuple[Clean, Clean, list[str]],
    settings: HybridSettings,
    known: KnownEventSites,
) -> dict[Site, dict[str, Any]]:
    """Per run site, its best within-run known-event minority (see module docstring)."""
    tally = _InRun(*observed)
    out: dict[Site, dict[str, Any]] = {}
    for site, per in tally.lengths.items():
        base, major = peers.meta[site]
        c: Counter[int] = Counter()
        for counter in per.values():
            c.update(counter)
        n = sum(tally.bounded[site].values())
        for minor, n_minor in c.items():
            if minor == major or n_minor < settings.phase_min_minor_reads:
                continue
            if (base, major, minor) not in known.runs:
                continue
            obs = _in_run_obs(peers, tally, observed, site, (major, minor), settings)
            share, bound = _bound(obs, settings)
            best = out.get(site)
            if bound >= settings.phase_run_minor_min_share and (
                best is None or bound > best["share_lower_bound"]
            ):
                out[site] = _entry(site, major, minor, (n_minor, n), (share, bound), IN_RUN)
    return out


def _entry(
    site: Site,
    major: int,
    minor: int,
    counts: tuple[int, int],
    bound: tuple[float, float],
    observation: str,
) -> dict[str, Any]:
    """A tier site record (``counts``: minor reads, reads; ``bound``: share, its bound)."""
    return {
        "site": site,
        "major": major,
        "minor": minor,
        "af": round(counts[0] / counts[1], AF_DECIMALS),
        "n": counts[1],
        "share": round(bound[0], AF_DECIMALS),
        "share_lower_bound": round(bound[1], AF_DECIMALS),
        "observation": observation,
    }


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
    *,
    impure: Clean | None = None,
) -> list[dict[str, Any]]:
    """Run sites whose minority length is significantly above the expected stutter.

    ``clean`` holds the peak's per-read clean run observations (``features(clean=...)``)
    and ``meta`` the site table's runs. For each run and each other observed length
    with at least ``phase_min_minor_reads`` reads within ``phase_run_bg_window`` bases
    of the modal length, the lower confidence bound of the deconvolved share must
    reach ``phase_run_minor_min_share``; with ``phase_run_minor_scope =
    "known_events"`` only lengths that make the change a dictionary template's
    site-table signature (``known_events``) are tested. With ``impure`` (the reads'
    bounded but impure run lengths, aligned with ``clean``) and
    ``phase_run_minor_in_run``, each run's known-event signature lengths are also tested
    inside the run (module docstring). Largest bound first; each run reports its best
    length and the ``observation`` it was found in (``CLEAN`` or ``IN_RUN``).
    """
    if settings.phase_run_minor_scope == "off":
        return []
    modal = clean_meta(clean, meta)
    peers = _Peers(clean, strands, modal)
    found: dict[Site, dict[str, Any]] = {}
    if impure is not None and settings.phase_run_minor_in_run:
        found = _in_run_sites(peers, (clean, impure, strands), settings, known)
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
        if best is None or best[1] <= found.get(site, {}).get("share_lower_bound", 0.0):
            continue
        counts = (c[best[2]], sum(c.values()))
        found[site] = _entry(site, major, best[2], counts, best[:2], CLEAN)
    return sorted(found.values(), key=lambda site: -site["share_lower_bound"])
