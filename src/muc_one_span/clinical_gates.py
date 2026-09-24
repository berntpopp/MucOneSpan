"""Evidence gates applied before a clinical decision banner is chosen.

Pure functions over summary dictionaries: they never read files or run tools.
Missing evidence is never treated as support.
"""

from __future__ import annotations

from typing import Any

SUPPORTED_VCF_STATUSES = frozenset({"exact_sequence_concordance"})
# Per-allele depth statuses that block a negative call and a PATHOGENIC carrier.
LOW_DEPTH_STATUSES = frozenset({"low", "insufficient"})
# read_support.status values a producer may emit; only "supported" is support.
READ_SUPPORT_STATUSES = frozenset(
    {"supported", "insufficient_depth", "discordant", "not_supported", "not_localized"}
)
_DEPTH_BASIS_LABELS = {
    "primary_alignment_records": "primary alignments",
    "spanning_reads": "spanning reads",
}


def mutation_supported(mutation: dict[str, Any]) -> bool:
    """Return True only for explicit sequence-level support (VCF or read evidence).

    ``read_support.status == "supported"`` is accepted for read-level evidence.
    A legacy ``vcf_support=True`` without a status keeps its historical meaning;
    a missing ``vcf_support`` key is not support.
    """
    read_support = mutation.get("read_support")
    if isinstance(read_support, dict) and read_support.get("status") == "supported":
        return True
    if mutation.get("vcf_support") is True:
        status = mutation.get("vcf_support_status")
        return status is None or status in SUPPORTED_VCF_STATUSES
    return False


LEGACY_MIN_TOTAL_READS = 30

_GENOTYPE_REASONS = {
    "heterozygous_within_length_partition": (
        "heterozygous call left within the length-partitioned allele; "
        "consensus uses unresolved (IUPAC) selection"
    ),
    "unresolved_genotype_records": "conflicting or incomplete genotype records",
    "residual_heterogeneity": (
        "residual read heterogeneity on this allele "
        "(possible unresolved mixture, chimera or mosaicism)"
    ),
}


def mutation_blockers(mutation: dict[str, Any]) -> list[str]:
    """List every reason an observed mutation cannot support a PATHOGENIC banner."""
    blockers: list[str] = []
    if mutation.get("frameshift") is not True:
        blockers.append("frameshift not established")
    if mutation.get("template_match") is not True or not mutation.get("mutation_name"):
        blockers.append("event identity not established (no exact dictionary template)")
    if mutation.get("localization_status") == "ambiguous":
        blockers.append("localization ambiguous")
    if not mutation_supported(mutation):
        read_support = mutation.get("read_support")
        status = mutation.get("vcf_support_status")
        if isinstance(read_support, dict):
            state = read_support.get("status") or "status unavailable"
            blockers.append(f"read-level support {state}")
        elif status == "heterozygous_genotype_unresolved":
            blockers.append("heterozygous genotype not resolved to one allele")
        else:
            blockers.append(
                f"no explicit sequence-level support ({status or 'status unavailable'})"
            )
    return blockers


def allele_gate_reasons(info: Any, label: str) -> list[str]:
    """Reasons an allele's selection, genotype, length or depth prevents a negative call."""
    if not isinstance(info, dict) or not info:
        return []
    reasons: list[str] = []
    selection = info.get("selection_status")
    if isinstance(selection, str) and selection.startswith("unresolved"):
        fraction = info.get("secondary_mode_fraction")
        detail = "" if fraction is None else f"; secondary mode fraction {fraction}"
        extra = f" {info['selection_detail']}" if info.get("selection_detail") else ""
        reasons.append(f"{label}: allele selection unresolved ({selection}{detail}).{extra}")
    length, reference_length = info.get("length"), info.get("reference_length")
    length_status = info.get("length_status")
    if length_status is None:  # Legacy summary without selection_qc: compare conservatively.
        differs = (
            isinstance(length, int) and isinstance(reference_length, int)
        ) and length != reference_length
    else:
        differs = length_status == "cluster_center_differs_from_consensus_contig"
    if differs:
        reasons.append(
            f"{label}: reported length {length} differs from the consensus contig length "
            f"{reference_length}."
        )
    if info.get("depth_status") in LOW_DEPTH_STATUSES:
        basis = info.get("depth_basis") or "primary_alignment_records"
        reasons.append(
            f"{label}: {info.get(basis)} {_DEPTH_BASIS_LABELS.get(basis, basis)}, below the "
            f"per-allele depth gate ({info.get('depth_threshold')})."
        )
    genotype = info.get("allele_genotype_status")
    if genotype in _GENOTYPE_REASONS:
        reasons.append(f"{label}: {_GENOTYPE_REASONS[genotype]}.")
    return reasons
