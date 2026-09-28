"""Split a length peak on its only heterozygous event.

``phase.split_by_linked_sites`` leaves a peak unsplit when fewer than
``min_linked_sites`` linked events support two haplotypes; with exactly one candidate
event that is ``unconfirmed_single_site`` and blocks a negative call. Two alleles of
the same repeat-unit length that differ at a single frameshift-sized event (a
one-base duplication at the end of a unit, say) have exactly one such event, so the
event itself must separate the haplotypes or the carrier allele is merged into the
wild-type consensus.

When ``phase_single_event_split`` allows it ("indel": the event changes the sequence
length, which every frameshift does; "all": any event), the members are split by
their allele at the event's most balanced site. A run-length site assigns each read to
the more likely of the two run lengths under the strand's stutter profile
(``run_strand``), so a stuttered observation is not simply dropped; a column or
insertion site uses the observed allele. Reads with neither allele (or a tie) are
returned in ``unassigned`` for the caller to reassign. The candidate already passed
``het_af_min``, the run-background floor and the strand-bias test; the split is also
refused when either group is below ``het_min_group`` of the members. The engine tries
this only for a length model with a single peak (``engine._groups``); with two length
peaks each peak already is one allele.

The split reads are selected by the event itself, so the allele's read support
(``evidence``) is conditional on the split. The split is therefore made only when a
peak-level gate independent of the split passes: the one-sided lower confidence bounds
(``phase_single_event_alpha``) of the stutter-deconvolved minor share **and** of the
major share, over a fresh seeded sample of at most ``phase_single_event_bound_reads``
reads, must both reach ``phase_single_event_min_share`` (at least ``het_af_min``; above
the share of a site-specific wild-type artefact, which the split cannot tell from a
real minor allele). Both sides, because the minor allele is the site's non-draft allele:
a draft that follows an artefact (a normal with a delinsAT-shaped artefact in 40% of
reads was drafted with it) makes the wild type the "minor" allele, and bounding only
that side split off the artefact group. For a run site the bounds must also reach the
floor with each length's own run-length stutter profile
(``run_minor.length_aware_share_bound``), besides the per-strand error profiles of
``run_strand``. Otherwise the peak stays ``unconfirmed_single_site`` (INCONCLUSIVE,
located).

Every tunable is a validated ``HybridSettings`` field taken from ``settings``.
"""

from __future__ import annotations

import random
from typing import Any

from muc_one_span.hybrid.phase import PhaseResult
from muc_one_span.hybrid.phase_sites import (
    GAP,
    Meta,
    Site,
    events,
    features,
    modal_meta,
    site_counts,
    top_site,
)
from muc_one_span.hybrid.run_minor import length_aware_share_bound
from muc_one_span.hybrid.run_strand import (
    ShiftProfile,
    length_prob,
    share_lower_bound,
    shift_profiles,
)
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

UNCONFIRMED = "unconfirmed_single_site"
SINGLE_EVENT = "single_event"


def is_indel(site: dict[str, Any]) -> bool:
    """True when the site's two alleles differ in length (a run length, gap or insert)."""
    kind = site["site"][0]
    if kind == "run":
        return bool(site["major"] != site["minor"])
    if kind == "ins":
        return len(site["major"]) != len(site["minor"])
    return GAP in (site["major"], site["minor"])


def _run_classifier(
    feats: list[dict[Site, Any]],
    strands: list[str],
    meta: Meta,
    site: dict[str, Any],
    s: HybridSettings,
) -> tuple[ShiftProfile, ShiftProfile]:
    key = site["site"]
    return (
        shift_profiles(feats, strands, meta, key, site["major"], s),
        shift_profiles(feats, strands, meta, key, site["minor"], s),
    )


def _allele(
    f: dict[Site, Any],
    strand: str,
    site: dict[str, Any],
    profiles: tuple[ShiftProfile, ShiftProfile] | None,
) -> int | None:
    """1 for the minor allele, 0 for the major, None when uninformative or tied."""
    observed = f.get(site["site"])
    if observed is None:
        return None
    if profiles is None:
        return 1 if observed == site["minor"] else 0 if observed == site["major"] else None
    p_major = length_prob(profiles[0], strand, observed, site["major"])
    p_minor = length_prob(profiles[1], strand, observed, site["minor"])
    return None if p_major == p_minor else int(p_minor > p_major)


def _bound_sample(n: int, s: HybridSettings) -> list[int]:
    """Indices of the share bound's fresh seeded sample (at most bound_reads reads)."""
    idx = list(range(n))
    if n > s.phase_single_event_bound_reads:
        idx = sorted(random.Random(s.seed).sample(idx, s.phase_single_event_bound_reads))
    return idx


def _length_aware_bounds(
    feats: list[dict[Site, Any]],
    strands: list[str],
    meta: Meta,
    site: dict[str, Any],
    s: HybridSettings,
) -> list[float]:
    """The share bounds of both run lengths with each length's own stutter profile.

    The run-minority tier's model (``run_minor.length_aware_share_bound``), for the
    minor share and for the major share; a side whose length is too rarely observed
    to be bounded is left out.
    """
    idx = _bound_sample(len(feats), s)
    runs = [{k: v for k, v in feats[i].items() if k[0] == site["site"][0]} for i in idx]
    sides = ((site["major"], site["minor"]), (site["minor"], site["major"]))
    bounds = [
        length_aware_share_bound(
            runs,
            [strands[i] for i in idx],
            meta,
            site["site"],
            alleles,
            s.phase_single_event_alpha,
            s,
        )
        for alleles in sides
    ]
    return [b for b in bounds if b is not None]


def _share_bounds(
    feats: list[dict[Site, Any]],
    strands: list[str],
    site: dict[str, Any],
    profiles: tuple[ShiftProfile, ShiftProfile] | None,
    s: HybridSettings,
) -> tuple[float, float]:
    """Lower confidence bounds of the minor share and of the major share.

    At most ``phase_single_event_bound_reads`` reads (a fresh sample drawn with
    ``random.Random(s.seed)``), so the bound's power does not grow with depth. Run
    reads contribute their likelihood under each run length (stutter profiles);
    column reads their allele (reads with neither allele are skipped). The minor
    allele is the site's non-draft allele; a draft that follows an artefact makes
    the artefact the major allele, so both shares are bounded and both must reach the
    floor: the smaller group of a split is then at least the floor, whichever allele
    the draft took.
    """
    obs: list[tuple[float, float]] = []
    for i in _bound_sample(len(feats), s):
        observed = feats[i].get(site["site"])
        if observed is None:
            continue
        if profiles is not None:
            obs.append(
                (
                    length_prob(profiles[1], strands[i], observed, site["minor"]),
                    length_prob(profiles[0], strands[i], observed, site["major"]),
                )
            )
        elif observed in (site["minor"], site["major"]):
            obs.append((float(observed == site["minor"]), float(observed == site["major"])))
    alpha = s.phase_single_event_alpha
    return share_lower_bound(obs, alpha), share_lower_bound([(b, a) for a, b in obs], alpha)


def split_single_event(
    cons: str, members: list[SpanRead], res: PhaseResult, settings: HybridSettings
) -> PhaseResult | None:
    """Two groups split on the peak's only candidate event, or None to keep it unsplit."""
    mode = settings.phase_single_event_split
    if mode == "off" or res.basis != UNCONFIRMED or not res.sites:
        return None
    feats, meta = features(cons, [m.seq for m in members], settings)
    if len(events(res.sites, meta)) != 1:
        return None
    if mode == "indel" and not any(is_indel(site) for site in res.sites):
        return None
    site = top_site(res.sites)
    strands = [m.strand for m in members]
    modal = modal_meta(site_counts(feats), meta)
    profiles = (
        _run_classifier(feats, strands, modal, site, settings) if site["site"][0] == "run" else None
    )
    floor = settings.phase_single_event_min_share
    if min(_share_bounds(feats, strands, site, profiles, settings)) < floor:
        return None
    if profiles is not None and any(
        # A run site's shares must also clear the floor with length-aware stutter
        # profiles: the profiles above pool every run of the base when the major
        # length has no peer, which under length-dependent stutter under-estimates
        # that run's own stutter and inflates the share.
        bound < floor
        for bound in _length_aware_bounds(feats, strands, meta, site, settings)
    ):
        return None
    groups: list[list[SpanRead]] = [[], []]
    unassigned: list[SpanRead] = []
    for m, f in zip(members, feats, strict=True):
        allele = _allele(f, m.strand, site, profiles)
        (unassigned if allele is None else groups[allele]).append(m)
    if min(len(g) for g in groups) < settings.het_min_group * len(members):
        return None
    return PhaseResult(groups, SINGLE_EVENT, res.sites, candidate=site, unassigned=unassigned)
