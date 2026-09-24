"""Evidence gates applied before a clinical decision banner is chosen.

Pure functions over summary dictionaries: they never read files or run tools.
Missing evidence is never treated as support.
"""

from __future__ import annotations

from typing import Any

SUPPORTED_VCF_STATUSES = frozenset({"exact_sequence_concordance"})


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
        status = mutation.get("vcf_support_status")
        blockers.append(
            "heterozygous genotype not resolved to one allele"
            if status == "heterozygous_genotype_unresolved"
            else f"no explicit sequence-level support ({status or 'status unavailable'})"
        )
    return blockers


def allele_gate_reasons(info: Any, label: str) -> list[str]:
    """Reasons an allele's selection, genotype, length or depth prevents a negative call."""
    if not isinstance(info, dict) or not info:
        return []
    reasons: list[str] = []
    selection = info.get("selection_status")
    if isinstance(selection, str) and selection.startswith("unresolved"):
        reasons.append(
            f"{label}: allele selection unresolved ({selection}; secondary mode fraction "
            f"{info.get('secondary_mode_fraction')})."
        )
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
    if info.get("depth_status") == "low":
        reasons.append(
            f"{label}: {info.get('primary_alignment_records')} primary alignments, below the "
            f"per-allele depth gate ({info.get('depth_threshold')})."
        )
    genotype = info.get("allele_genotype_status")
    if genotype in _GENOTYPE_REASONS:
        reasons.append(f"{label}: {_GENOTYPE_REASONS[genotype]}.")
    return reasons
