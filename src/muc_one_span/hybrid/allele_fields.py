"""Per-allele output contract and sample-level selection status of the hybrid engine.

Every allele record keeps the ladder's allele contract (``length``, ``reads``,
``depth_status`` ...) so the shared gates in ``clinical_gates``/``report`` apply
unchanged, and adds hybrid provenance fields. Any status other than ``resolved``
starts with ``unresolved`` and therefore blocks a negative call
(``clinical_gates.allele_gate_reasons``).
"""

from __future__ import annotations

from typing import Any

from muc_one_span.hybrid.lengths import LengthModel
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

# A sample carries at most two MUC1 alleles (diploid autosomal locus): a structural
# domain constant, not a tunable.
PLOIDY = 2
# Split basis -> phase_status. "length" and "linked_sites" are independent haplotype
# evidence; an unconfirmed split leaves the peak unsplit and its status unresolved.
PHASE_STATUS = {
    "length": "phased",
    "linked_sites": "phased",
    "none": "no_informative_heterozygosity",
    "unconfirmed_single_site": "unresolved_single_site",
    "unconfirmed_group_size": "unresolved_group_size",
    "single_event": "phased_single_event",
    "unconfirmed_run_site": "unresolved_run_site",
}
# A single-event split (Task 15e) separates the reads by their allele at the only
# heterozygous event, so it is read-level haplotype evidence like a linked-site split.
INDEPENDENT_BASES = frozenset({"length", "linked_sites", "single_event"})
# Sample-level selection statuses of the hybrid engine (all but RESOLVED gate).
RESOLVED = "resolved"
MAX_ALLELES = "unresolved_max_alleles"
SINGLE_SITE = "unresolved_single_site"
GROUP_SIZE = "unresolved_group_size"
# Task 15g: an unsplit equal-length peak with a run above the lower safety floor.
RUN_SITE = "unresolved_run_site"
REJECTED_PEAK = "unresolved_rejected_peak"
UNASSIGNED_SPANNING = "unresolved_unassigned_spanning"
# Length-model rejection reason that means "a third allele-like peak" (lengths.py).
MAX_ALLELES_REASON = "max_alleles"
UNCONFIRMED_SPLIT_STATUS = {
    "unconfirmed_single_site": SINGLE_SITE,
    "unconfirmed_group_size": GROUP_SIZE,
    "unconfirmed_run_site": RUN_SITE,
}


def depth_status(spanning: int, h: HybridSettings) -> str:
    """adequate / low / insufficient from the allele's spanning-read count."""
    if spanning >= h.depth_adequate_spanning:
        return "adequate"
    return "low" if spanning >= h.depth_low_spanning else "insufficient"


def selection_status(
    model: LengthModel,
    split_bases: list[str],
    n_groups: int,
    unassigned_fraction: float,
    h: HybridSettings,
) -> str:
    """Sample-level selection status; any unresolved state blocks a negative call.

    Precedence: more than two allele groups (or a rejected third peak), an unconfirmed
    linked-site split, an unsplit peak with a run above the safety floor
    (``unconfirmed_run_site``), any other gate-relevant rejected length peak (including
    ``smear_ambiguous``), then too many spanning reads assigned to no allele.
    """
    reasons = {r["reason"] for r in model.gate_relevant_rejections}
    if n_groups > PLOIDY or MAX_ALLELES_REASON in reasons:
        return MAX_ALLELES
    for basis, status in UNCONFIRMED_SPLIT_STATUS.items():
        if basis in split_bases:
            return status
    if reasons:
        return REJECTED_PEAK
    if unassigned_fraction > h.max_unassigned_spanning_fraction:
        return UNASSIGNED_SPANNING
    return RESOLVED


def selection_detail(
    model: LengthModel,
    n_groups: int,
    unassigned: int,
    unassigned_fraction: float,
    *,
    unresolved_sites: list[str],
) -> str:
    """Human-readable reason for the selection status (shown next to the gate reason).

    ``unresolved_sites`` locates each heterozygous site that left a peak unsplit
    (``engine.located_site``), so the gate reason names the repeat to inspect.
    """
    rejected = ", ".join(
        f"{r['units']} units ({r['reason']}, {r['support']} reads)"
        for r in model.gate_relevant_rejections
    )
    return (
        f"{n_groups} allele group(s); gate-relevant rejected length peaks: "
        f"{rejected or 'none'}; {unassigned} spanning read(s) assigned to no allele "
        f"({unassigned_fraction:.1%})."
        + (f" Unresolved: {'; '.join(unresolved_sites)}." if unresolved_sites else "")
    )


def allele_info(
    name: str,
    seq: str,
    members: list[SpanRead],
    assigned: int,
    basis: str,
    residual: list[dict[str, Any]],
    status: tuple[str, str],
    unit_bp: int,
    fixed_repeat_count: int,
    h: HybridSettings,
    *,
    concordance: float,
) -> dict[str, Any]:
    """One allele's record; ``status`` is the sample (selection_status, selection_detail).

    ``length`` counts every repeat unit of the motif 1..9 consensus (as the ladder's
    ``length``); ``canonical_repeats`` excludes the reference layout's fixed repeats.

    ``concordance`` (``hybrid.polish.consensus_concordance`` on the final polished
    consensus and the reads that built it) is the hybrid engine's own read-support
    evidence and is reported as ``consensus_concordance_fraction``, alongside a
    ``classification_confidence_status`` marker: the ladder's ``classify.py``
    ``confidence``/``allele_confidence`` (``repeats.json``, the CLI's printed
    ``confidence:`` line, and the report's "Allele confidence" tile) is a
    dictionary-fit heuristic computed identically for both engines and carries no
    hybrid reconstruction evidence, so it does not apply here.
    """
    units = round(len(seq) / unit_bp)
    depth = depth_status(len(members), h)
    selection, detail = status
    return {
        "engine": "hybrid",
        "length": units,
        "canonical_repeats": units - fixed_repeat_count,
        "fixed_repeat_count": fixed_repeat_count,
        "reference_length": units,
        "length_status": "consistent_with_consensus_contig",
        "length_basis": "spanning_peak",
        "consensus_basis": "poa",
        "reads": len(members),
        "spanning_reads": len(members),
        "assigned_reads": assigned,
        "depth_status": depth,
        "depth_basis": "spanning_reads",
        "depth_threshold": h.depth_adequate_spanning,
        "selection_status": selection,
        "secondary_mode_fraction": None,
        "selection_detail": detail,
        "split_basis": basis,
        "phase_status": PHASE_STATUS[basis],
        "independent_haplotype_evidence": basis in INDEPENDENT_BASES,
        "residual_sites": residual,
        "heterozygous_sites": residual,
        "variant_filter": None,
        "allele_genotype_status": (
            "residual_heterogeneity" if residual else "not_applicable_read_consensus"
        ),
        "reconstruction_status": (
            "complete_segmentation" if depth == "adequate" else "read_consensus_low_depth"
        ),
        "sequence_source": f"hybrid:{name}",
        "contig_name": f"hybrid_{name}",
        "vcf_path": None,
        "consensus_concordance_fraction": concordance,
        "classification_confidence_status": "not_applicable_dictionary_fit_heuristic",
    }
