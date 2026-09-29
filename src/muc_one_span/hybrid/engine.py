"""Hybrid engine orchestration: spans -> lengths -> POA -> phase -> assign -> polish.

Returns the existing allele contract (plus added fields) and a separate ``block`` that
the pipeline stores as ``summary["hybrid"]``; nothing engine-specific is stored inside
``alleles`` except per-allele fields. Deviations from spec §3 in this version: no
ladder-assisted length prior or ladder-seeded consensus (S2/S3), depth is judged on
spanning reads only, and read-level support uses the spanning members only.

Every tunable is read from the ``RuntimeSettings`` passed in (``settings.hybrid``) and
handed to each stage explicitly; the repeat-unit length comes from the repeat
dictionary and the fixed-repeat count from ``settings.reference_layout``.
"""

from __future__ import annotations

import json
import random
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TypeVar

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.allele_fields import (
    PLOIDY,
    RESOLVED,
    SINGLE_SITE,
    UNCONFIRMED_SPLIT_STATUS,
    allele_info,
    selection_detail,
    selection_status,
    single_group_fields,
)
from muc_one_span.hybrid.assign import (
    OFF_TARGET,
    UNDECIDED,
    assign_read,
    assign_reads,
    hybrid_references,
    trim_to_draft,
)
from muc_one_span.hybrid.evidence import (
    FRACTION_DECIMALS,
    event_read_support,
    residual_sites,
)
from muc_one_span.hybrid.known_events import (
    KnownEventSites,
    known_event_sites,
    warn_blind_templates,
)
from muc_one_span.hybrid.lengths import LengthModel, fit_length_model
from muc_one_span.hybrid.phase import PhaseResult, split_by_linked_sites
from muc_one_span.hybrid.phase_groups import explained_group
from muc_one_span.hybrid.poa import PoaBackend, get_backend
from muc_one_span.hybrid.polish import consensus_concordance, draft_consensus, polish
from muc_one_span.hybrid.reads_io import extra_versions, read_input
from muc_one_span.hybrid.single_event import SINGLE_EVENT, split_single_event
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.run_status import InsufficientEvidenceError
from muc_one_span.settings import HybridSettings, RuntimeSettings

__all__ = [
    "HybridResult",
    "annotate_read_support",
    "extra_versions",
    "read_input",
    "reconstruct_alleles",
]
T = TypeVar("T")
# Split bases of an unsplit equal-length peak held back by a safety tier.
RUN_SITE_BASIS = "unconfirmed_run_site"
BIASED_SITE_BASIS = "unconfirmed_strand_biased_site"
# Split basis of an unsplit peak (any peak count) held back by the Task 15l run-minority
# tier.
RUN_MINOR_BASIS = "unconfirmed_run_minor"


@dataclass
class HybridResult:
    """Engine output; ``block`` becomes ``summary["hybrid"]`` and is never put in alleles."""

    alleles: dict[str, Any]
    consensus_paths: dict[str, Path]
    block: dict[str, Any]
    members: dict[str, list[tuple[str, str]]]


@dataclass
class _Group:
    """Spanning members of one allele candidate, its split basis and its POA draft."""

    members: list[SpanRead]
    basis: str
    draft: str


@dataclass
class _Phased:
    """The allele groups of every length peak and what the phase stage recorded.

    ``split_bases`` holds each peak's split basis, ``leftover`` the spanning reads no
    group took, ``located`` the located sites (``located_site``) keyed by split basis,
    ``dropped`` the sites dropped as explained by a low-accuracy read subset
    (``quality_site``), ``quality_split`` those among them whose removal let an
    equal-length peak be split on its remaining event, and ``excluded`` the linked
    groups of a two-peak model not counted as alleles (``phase_groups``; both Task 15k).
    """

    groups: list[_Group] = field(default_factory=list)
    split_bases: list[str] = field(default_factory=list)
    leftover: list[SpanRead] = field(default_factory=list)
    located: dict[str, list[str]] = field(default_factory=dict)
    dropped: list[dict[str, Any]] = field(default_factory=list)
    quality_split: list[dict[str, Any]] = field(default_factory=list)
    excluded: list[dict[str, Any]] = field(default_factory=list)


def _cap(items: list[T], limit: int, rng: random.Random) -> list[T]:
    """At most ``limit`` items, sampled with the seeded RNG (input order is arbitrary)."""
    return items if len(items) <= limit else rng.sample(items, limit)


def _draft(members: list[SpanRead], h: HybridSettings, rng: random.Random, be: PoaBackend) -> str:
    return draft_consensus(
        members,
        h.n_poa,
        rng,
        be,
        sample_window_floor_bp=h.poa_sample_window_floor_bp,
        sample_window_frac=h.poa_sample_window_frac,
    )


def _reassign(
    unassigned: list[SpanRead], groups: list[_Group], h: HybridSettings
) -> list[SpanRead]:
    """Place phase-unassigned spanning reads by edit distance to the group drafts (S6).

    A read joins a group only when that draft wins by ``assign_margin`` in the read's
    own (forward) orientation; everything else is returned and counted as a spanning
    read assigned to no allele, never dropped.
    """
    refs = {str(i): g.draft for i, g in enumerate(groups)}
    leftover = []
    for sp in unassigned:
        got = assign_read(sp.seq, refs, h.assign_margin, h.assign_max_error_rate)
        if got.allele in refs and got.oriented == sp.seq:
            groups[int(got.allele)].members.append(sp)
        else:
            leftover.append(sp)
    return leftover


def _single_event(
    draft: str,
    members: list[SpanRead],
    split: PhaseResult,
    h: HybridSettings,
    rng: random.Random,
    backend: PoaBackend,
) -> tuple[PhaseResult, list[_Group] | None]:
    """Promote an unconfirmed single-event peak to a split when both drafts differ."""
    promoted = split_single_event(draft, members, split, h)
    if promoted is None:
        return split, None
    sub = [_Group(g, promoted.basis, _draft(g, h, rng, backend)) for g in promoted.groups]
    if sub[0].draft == sub[1].draft:
        return split, None
    return promoted, sub


def _run_site_tier(split: PhaseResult, single_peak: bool) -> PhaseResult:
    """Keep an unsplit peak from a negative call on a sub-floor site.

    For a single-peak model, a peak with no candidate site ("none") whose runs include
    one above the lower safety floor (``PhaseResult.run_excess``) becomes
    ``unconfirmed_run_site``; one with a column or insertion site refused only for
    strand bias (``PhaseResult.strand_biased``, Task 15i) becomes
    ``unconfirmed_strand_biased_site``. In any peak (Task 15l), a run whose minority
    length is significantly above its expected stutter (``PhaseResult.run_minor``)
    makes it ``unconfirmed_run_minor``. Each stays one group with no event, but its
    selection status blocks a negative call and the site (largest excess, allele
    fraction or share bound) is named in the reason.
    """
    if split.basis != "none":
        return split
    tiers = [(RUN_SITE_BASIS, split.run_excess), (BIASED_SITE_BASIS, split.strand_biased)]
    for basis, sites in [*(tiers if single_peak else []), (RUN_MINOR_BASIS, split.run_minor)]:
        if sites:
            return PhaseResult(split.groups, basis, sites, candidate=sites[0])
    return split


def located_site(candidate: dict[str, Any], unit_bp: int) -> str:
    """Gate reason naming an unresolved heterozygous site by its 1-based repeat unit."""
    kind, pos = candidate["site"]
    return (
        f"unresolved heterozygous site at repeat {pos // unit_bp + 1} "
        f"({kind} {candidate['major']!r}>{candidate['minor']!r}, AF {candidate['af']})"
    )


def quality_site(site: dict[str, Any], unit_bp: int) -> dict[str, Any]:
    """Report record of a site dropped as explained by a low-accuracy read subset."""
    kind, pos = site["site"]
    return {
        "repeat": pos // unit_bp + 1,
        "kind": kind,
        **{
            k: site[k]
            for k in ("major", "minor", "af", "af_high_quality", "quality_p", "af_bound_p")
        },
    }


def _groups(
    model: LengthModel,
    h: HybridSettings,
    rng: random.Random,
    backend: PoaBackend,
    unit_bp: int,
    insertions: KnownEventSites,
) -> _Phased:
    """Split each length peak by linked sites (or its single event).

    The low-accuracy-subset rule (Task 15j) applies only when the length model found
    ``PLOIDY`` peaks: each allele then has its own peak, and a within-peak site is not
    a further allele.

    A single-event split is tried only when the length model found fewer than
    ``PLOIDY`` peaks (the equal-length heterozygote): with two length peaks each peak
    already is one allele, so a within-peak single-site mixture is not a further
    haplotype, and splitting it would only move stutter or error reads out of an allele
    and make that allele's read support circular. Such a peak stays unconfirmed. The
    15g run-site and 15i strand-biased tiers (``_run_site_tier``) apply under the same
    condition; its Task 15l run-minority tier applies to every unsplit peak. In a
    two-peak model a peak held back only by that tier keeps the split basis
    ``length`` for its allele (the allele is still phased by its length peak); the
    sample's selection status, which every allele carries and which gates a negative
    call, is ``unresolved_run_minor`` through ``split_bases``. When
    more than one event leaves such a peak unconfirmed and the sites the 15j test keeps
    form exactly one event (``phase_quality_single_event``, Task 15k), the split is
    tried on that event under the same gates; the dropped sites are then recorded in
    ``quality_split`` and the caller keeps the selection unresolved. Neither quality
    rule drops a site carrying a dictionary event signature (``insertions``).
    """
    out = _Phased()
    two_peaks = len(model.peaks) == PLOIDY
    for peak in model.peaks:
        draft = _draft(peak.members, h, rng, backend)
        split = split_by_linked_sites(
            draft,
            peak.members,
            h,
            rng,
            quality_filter=two_peaks,
            single_event_quality=not two_peaks and h.phase_quality_single_event,
            insertions=insertions,
        )
        sub: list[_Group] | None = None
        if len(model.peaks) < PLOIDY:
            alternative = split.quality_single_event
            split, sub = _single_event(draft, peak.members, split, h, rng, backend)
            if sub is None and alternative is not None:
                promoted, sub = _single_event(draft, peak.members, alternative, h, rng, backend)
                if sub is not None:
                    split = promoted
                    split.quality_associated = alternative.quality_associated
                    out.quality_split.extend(
                        quality_site(site, unit_bp) for site in alternative.quality_associated
                    )
        split = _run_site_tier(split, single_peak=not two_peaks)
        out.dropped.extend(quality_site(site, unit_bp) for site in split.quality_associated)
        out.split_bases.append(split.basis)
        if split.candidate is not None and (sub is not None or len(split.groups) == 1):
            out.located.setdefault(split.basis, []).append(located_site(split.candidate, unit_bp))
        if sub is None and len(split.groups) == 1:
            by_length = two_peaks and split.basis in ("none", RUN_MINOR_BASIS)
            out.groups.append(
                _Group(split.groups[0], "length" if by_length else split.basis, draft)
            )
            continue
        if sub is None:
            sub = [_Group(g, split.basis, _draft(g, h, rng, backend)) for g in split.groups if g]
        out.leftover.extend(_reassign(split.unassigned, sub, h))
        if two_peaks and split.basis == "linked_sites":
            sub = _exclude_explained(sub, out, h, unit_bp)
        out.groups.extend(sub)
    return out


def _exclude_explained(
    sub: list[_Group], out: _Phased, h: HybridSettings, unit_bp: int
) -> list[_Group]:
    """Drop a linked group explained by low-accuracy reads from the allele count.

    Its reads join no allele and are counted with the reads assigned to no allele
    (``out.leftover``); the evidence goes to ``out.excluded`` (``phase_groups``).
    """
    hit = explained_group([g.members for g in sub], [g.draft for g in sub], h, unit_bp)
    if hit is None:
        return sub
    index, record = hit
    out.leftover.extend(sub[index].members)
    out.excluded.append(record)
    return [g for i, g in enumerate(sub) if i != index]


def quality_split_detail(sites: list[dict[str, Any]]) -> str:
    """Selection-detail sentence for a single-event split made after dropping sites."""
    named = "; ".join(
        f"repeat {s['repeat']} ({s['kind']} {s['major']!r}>{s['minor']!r}, AF {s['af']}, "
        f"high-quality AF {s['af_high_quality']})"
        for s in sites
    )
    return (
        " The equal-length peak was split on its single event after dropping site(s) "
        f"explained by low-accuracy reads: {named}; a negative call stays blocked."
    )


def _write_allele(output_dir: Path, name: str, cons: str, n_reads: int) -> Path:
    path = output_dir / f"consensus_{name}.fa"
    path.write_text(f">{name}\n{cons}\n")
    context = {"engine": "hybrid", "reads": n_reads, "full_consensus_path": str(path)}
    (output_dir / f"consensus_{name}_context.json").write_text(
        json.dumps({**context, "vcf_path": None}, indent=2) + "\n"
    )
    return path


def reconstruct_alleles(
    input_path: Path, output_dir: Path, rd: RepeatDictionary, settings: RuntimeSettings
) -> HybridResult:
    """Reconstruct up to two allele sequences with explicit evidence and uncertainty."""
    h = settings.hybrid
    rng = random.Random(h.seed)
    backend = get_backend(h.poa_backend)
    anchors = Anchors.from_dictionary(rd, h, settings.reference_layout)
    unit_bp = anchors.unit_bp
    cats = categorize_reads(read_input(input_path), anchors, h)
    model = fit_length_model(cats.spanning, h, anchors)
    if not model.peaks:
        raise InsufficientEvidenceError("hybrid: no allele length peak passed the thresholds")
    known = known_event_sites(rd, h)
    warn_blind_templates(known, h)
    phased = _groups(model, h, rng, backend, unit_bp, known)
    groups, split_bases, phase_leftover = phased.groups, phased.split_bases, phased.leftover
    located, quality_dropped = phased.located, phased.dropped
    unresolved = [site for basis in UNCONFIRMED_SPLIT_STATUS for site in located.get(basis, [])]
    n_unassigned = len(model.unassigned) + len(phase_leftover)
    unassigned_fraction = n_unassigned / model.total
    selection = selection_status(model, split_bases, len(groups), unassigned_fraction, h)
    detail = selection_detail(
        model, len(groups), n_unassigned, unassigned_fraction, unresolved_sites=unresolved
    )
    if phased.quality_split:
        # Task 15k: the split may turn INCONCLUSIVE into PATHOGENIC, never NEGATIVE.
        selection = SINGLE_SITE if selection == RESOLVED else selection
        detail += quality_split_detail(phased.quality_split)
    ranked = sorted(groups, key=lambda g: -len(g.members))
    kept, dropped = ranked[:PLOIDY], ranked[PLOIDY:]
    kept.sort(key=lambda g: statistics.median(m.length for m in g.members))  # allele_1 shorter
    names = [f"allele_{i + 1}" for i in range(len(kept))]
    refs = hybrid_references(
        {n: g.draft for n, g in zip(names, kept, strict=True)}, rd, h.assign_flank_bp
    )
    extra = assign_reads(
        cats.left_anchored + cats.right_anchored + cats.internal_or_offtarget, refs, h
    )
    alleles: dict[str, Any] = {}
    paths: dict[str, Path] = {}
    members: dict[str, list[tuple[str, str]]] = {}
    seqs: dict[str, str] = {}
    fixed = settings.reference_layout.fixed_repeat_count
    for name, group in zip(names, kept, strict=True):
        partial = [
            trim_to_draft(r, refs[name], h.assign_flank_bp, len(group.draft))
            for r in _cap(extra[name], h.polish_max_reads, rng)
        ]
        full_reads = [m.seq for m in _cap(group.members, h.polish_max_reads, rng)]
        partial_reads = [p for p in partial if len(p) >= h.polish_partial_min_units * unit_bp]
        cons, info = polish(
            group.draft,
            full_reads,
            partial_reads,
            rounds=h.polish_rounds,
            hp_vote=h.hp_vote,
            insertion_majority_frac=h.polish_insertion_majority_frac,
            hp_min_run=h.hp_vote_min_run,
        )
        # Real hybrid read-support evidence: reuses the exact reads that built/polished
        # ``cons`` (no extra rng.sample draws, which would disturb downstream
        # determinism) rather than the ladder's dictionary-fit classify.py confidence.
        concordance = consensus_concordance(cons, full_reads, partial_reads)
        qc_reads = [m.seq for m in _cap(group.members, h.qc_residual_max_reads, rng)]
        residual = residual_sites(cons, qc_reads, h.qc_residual_af, min_run=h.qc_residual_min_run)
        alleles[name] = allele_info(
            name,
            cons,
            group.members,
            len(extra[name]),
            group.basis,
            residual,
            (selection, detail),
            unit_bp,
            fixed,
            h,
            concordance=round(concordance, FRACTION_DECIMALS),
        )
        alleles[name]["polish"] = info
        members[name] = [(m.seq, m.strand) for m in group.members]
        paths[name] = _write_allele(output_dir, name, cons, len(group.members))
        seqs[name] = cons
    if SINGLE_EVENT in split_bases and len(set(seqs.values())) < len(seqs):
        # A single-event split whose polished alleles are identical did not separate
        # the haplotypes: keep the negative call blocked, as for an unsplit site.
        selection = SINGLE_SITE
        detail += (
            " The single-event split gave identical allele sequences. Unresolved: "
            f"{'; '.join(located[SINGLE_EVENT])}."
        )
        for n in names:
            alleles[n].update(selection_status=selection, selection_detail=detail)
    residual_any = any(alleles[n]["residual_sites"] for n in names)
    if len(kept) == 1:
        alleles.update(single_group_fields(alleles["allele_1"], selection, residual_any))
    else:
        lengths = [alleles[n]["length"] for n in names]
        alleles.update(
            homozygous=False,
            same_length=len(set(lengths)) == 1,
            allele_multiplicity_status="resolved",
        )
    alleles["observed_length_candidates"] = [round(p.center_bp / unit_bp) for p in model.peaks]
    # Keep the historical key order (multiplicity last).
    alleles["allele_multiplicity_status"] = alleles.pop("allele_multiplicity_status")
    block = {
        "engine": "hybrid",
        "assay": settings.run.assay,
        "poa_backend": backend.name,
        "unit_bp": unit_bp,
        "read_categories": cats.counts(),
        "rejected_peaks": model.rejected,
        "short_product_fraction": round(model.short_product_fraction, FRACTION_DECIMALS),
        "dimer_product_reads": len(model.dimer_products),
        "dimer_product_fraction": round(model.dimer_product_fraction, FRACTION_DECIMALS),
        "unassigned_spanning_reads": n_unassigned,
        "phase_unassigned_spanning_reads": len(phase_leftover),
        "unassigned_spanning_fraction": round(unassigned_fraction, FRACTION_DECIMALS),
        "undecided_reads": len(extra[UNDECIDED]),
        "off_target_reads": len(extra[OFF_TARGET]),
        "selection_status": selection,
        "selection_detail": detail,
        "split_bases": split_bases,
        "quality_associated_sites": quality_dropped,
        "quality_excluded_groups": phased.excluded,
        "dropped_groups": [
            {
                "spanning_reads": len(g.members),
                "median_units": round(statistics.median(m.length for m in g.members) / unit_bp),
                "split_basis": g.basis,
            }
            for g in dropped
        ],
    }
    (output_dir / "hybrid_reads.json").write_text(json.dumps(block, indent=2) + "\n")
    (output_dir / "hybrid_references.fa").write_text(
        "".join(f">hybrid_{k}\n{v}\n" for k, v in refs.items())
    )
    return HybridResult(alleles, paths, block, members)


def annotate_read_support(
    allele_key: str,
    sequence: str,
    result: dict[str, Any],
    *,
    rd: RepeatDictionary,
    members: dict[str, list[tuple[str, str]]],
    settings: HybridSettings,
) -> dict[str, Any]:
    """Attach read-level support; VCF support does not apply to a read consensus.

    ``members`` holds, per allele, the spanning reads assigned to that allele, already
    oriented to its consensus (the VNTR forward strand), with their original strand.
    """
    support = event_read_support(sequence, result, members.get(allele_key, []), rd, settings)
    for idx, mutation in enumerate(result.get("mutations_detected", [])):
        read_support = support.get(idx) or {"kind": "none", "status": "not_localized"}
        mutation["read_support"] = read_support
        mutation["vcf_support"] = False
        mutation["vcf_support_status"] = "not_applicable_read_consensus"
        mutation["support_status"] = f"read_level_{read_support['status']}"
    return result
