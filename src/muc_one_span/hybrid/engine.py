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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeVar

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.allele_fields import (
    PLOIDY,
    RESOLVED,
    allele_info,
    selection_detail,
    selection_status,
)
from muc_one_span.hybrid.assign import (
    OFF_TARGET,
    UNDECIDED,
    assign_read,
    assign_reads,
    hybrid_references,
    trim_to_draft,
)
from muc_one_span.hybrid.evidence import event_read_support, residual_sites
from muc_one_span.hybrid.lengths import LengthModel, fit_length_model
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.poa import PoaBackend, get_backend
from muc_one_span.hybrid.polish import consensus_concordance, draft_consensus, polish
from muc_one_span.hybrid.reads_io import extra_versions, read_input
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
# Output format precision of reported fractions (statuses use unrounded values).
FRACTION_DECIMALS = 4
T = TypeVar("T")


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


def _groups(
    model: LengthModel, h: HybridSettings, rng: random.Random, backend: PoaBackend
) -> tuple[list[_Group], list[str], list[SpanRead]]:
    """Split each length peak by linked sites; return groups, split bases, leftovers."""
    groups: list[_Group] = []
    split_bases: list[str] = []
    leftover: list[SpanRead] = []
    for peak in model.peaks:
        draft = _draft(peak.members, h, rng, backend)
        split = split_by_linked_sites(draft, peak.members, h, rng)
        split_bases.append(split.basis)
        if len(split.groups) == 1:
            two_peaks = split.basis == "none" and len(model.peaks) == PLOIDY
            groups.append(_Group(split.groups[0], "length" if two_peaks else split.basis, draft))
            continue
        sub = [_Group(g, split.basis, _draft(g, h, rng, backend)) for g in split.groups if g]
        leftover.extend(_reassign(split.unassigned, sub, h))
        groups.extend(sub)
    return groups, split_bases, leftover


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
    anchors = Anchors.from_dictionary(rd, h)
    unit_bp = anchors.unit_bp
    cats = categorize_reads(read_input(input_path), anchors, h)
    model = fit_length_model(cats.spanning, h, unit_bp)
    if not model.peaks:
        raise InsufficientEvidenceError("hybrid: no allele length peak passed the thresholds")
    groups, split_bases, phase_leftover = _groups(model, h, rng, backend)
    n_unassigned = len(model.unassigned) + len(phase_leftover)
    unassigned_fraction = n_unassigned / model.total
    selection = selection_status(model, split_bases, len(groups), unassigned_fraction, h)
    detail = selection_detail(model, len(groups), n_unassigned, unassigned_fraction)
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
    residual_any = any(alleles[n]["residual_sites"] for n in names)
    homozygous = len(kept) == 1 and selection == RESOLVED and not residual_any
    if len(kept) == 1:
        alleles["allele_2"] = {
            **alleles["allele_1"],
            "candidate_duplicate_of": "allele_1",
            "reconstruction_status": "not_separately_resolved",
            "independent_haplotype_evidence": homozygous,
        }
        alleles["sequence_identity_status"] = "resolved" if homozygous else "unresolved"
    lengths = [alleles[n]["length"] for n in names]
    alleles.update(
        {
            "homozygous": homozygous,
            "same_length": len(set(lengths)) == 1,
            "observed_length_candidates": [round(p.center_bp / unit_bp) for p in model.peaks],
            "allele_multiplicity_status": "resolved" if len(kept) == PLOIDY else "unresolved",
        }
    )
    block = {
        "engine": "hybrid",
        "assay": settings.run.assay,
        "poa_backend": backend.name,
        "unit_bp": unit_bp,
        "read_categories": cats.counts(),
        "rejected_peaks": model.rejected,
        "short_product_fraction": round(model.short_product_fraction, FRACTION_DECIMALS),
        "unassigned_spanning_reads": n_unassigned,
        "phase_unassigned_spanning_reads": len(phase_leftover),
        "unassigned_spanning_fraction": round(unassigned_fraction, FRACTION_DECIMALS),
        "undecided_reads": len(extra[UNDECIDED]),
        "off_target_reads": len(extra[OFF_TARGET]),
        "selection_status": selection,
        "selection_detail": detail,
        "split_bases": split_bases,
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
