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
