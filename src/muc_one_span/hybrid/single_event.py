"""Split a length peak on its only heterozygous event (Task 15e).

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
peak-level gate independent of the split passes: the one-sided lower confidence bound
(``phase_single_event_alpha``) of the stutter-deconvolved minor share, over a fresh
seeded sample of at most ``phase_single_event_bound_reads`` reads, must reach
``het_af_min``. Otherwise the peak stays ``unconfirmed_single_site`` (INCONCLUSIVE,
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


def _share_bound(
    feats: list[dict[Site, Any]],
    strands: list[str],
    site: dict[str, Any],
    profiles: tuple[ShiftProfile, ShiftProfile] | None,
    s: HybridSettings,
) -> float:
    """Lower confidence bound of the minor share over a fixed-size read sample.

    At most ``phase_single_event_bound_reads`` reads (a fresh sample drawn with
    ``random.Random(s.seed)``), so the bound's power does not grow with depth. Run reads contribute their
    likelihood under each run length (stutter profiles); column reads their allele
    (reads with neither allele are skipped).
    """
    idx = list(range(len(feats)))
    if len(idx) > s.phase_single_event_bound_reads:
        idx = sorted(random.Random(s.seed).sample(idx, s.phase_single_event_bound_reads))
    obs: list[tuple[float, float]] = []
    for i in idx:
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
    return share_lower_bound(obs, s.phase_single_event_alpha)


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
    if _share_bound(feats, strands, site, profiles, settings) < settings.het_af_min:
        return None
    groups: list[list[SpanRead]] = [[], []]
    unassigned: list[SpanRead] = []
    for m, f in zip(members, feats, strict=True):
        allele = _allele(f, m.strand, site, profiles)
        (unassigned if allele is None else groups[allele]).append(m)
    if min(len(g) for g in groups) < settings.het_min_group * len(members):
        return None
    return PhaseResult(groups, SINGLE_EVENT, res.sites, candidate=site, unassigned=unassigned)
