"""Evaluator-side causes for an incomplete reconstruction, kept per evaluated sample.

`reconstruction_flags` names why ``evaluate_sample`` could not credit a clean
diploid reconstruction, in a fixed order:

- the observation status when it is not ``completed`` (for example
  ``ambiguous_reconstruction``, ``insufficient_evidence``, ``execution_failed``);
- ``iupac_bases``: a predicted consensus holds a base outside ``ACGT`` (the
  same test ``artifacts.load_observation`` uses to mark the run ambiguous);
- ``unresolved_allele_alias``: an allele was reported only as a duplicate of
  another (``unresolved_allele_alias`` loader warning);
- ``producer_status_unresolved``: the run is ambiguous for neither reason
  above, so a producer phase or reconstruction status was not accepted;
- ``missing_allele``, ``unproven_duplicate_allele``, ``extra_allele``: the
  scoring row's allele cardinality counts are nonzero.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from .models import RunObservation

UNAMBIGUOUS_BASES = frozenset("ACGT")  # format definition: unambiguous nucleotide codes
_ALIAS_WARNING = "unresolved_allele_alias:"  # loader warning prefix (artifacts.load_observation)
_COUNT_FLAGS = (
    ("missing_alleles", "missing_allele"),
    ("unproven_duplicate_alleles", "unproven_duplicate_allele"),
    ("extra_alleles", "extra_allele"),
)


def reconstruction_flags(observation: RunObservation, row: Mapping[str, Any]) -> list[str]:
    """Ordered evaluator causes for ``observation`` and its scoring ``row`` (see module doc)."""
    flags: list[str] = []
    if observation.status != "completed":
        flags.append(observation.status)
    iupac = any(set(p.sequence) - UNAMBIGUOUS_BASES for p in observation.predictions)
    alias = any(w.startswith(_ALIAS_WARNING) for w in observation.warnings)
    if iupac:
        flags.append("iupac_bases")
    if alias:
        flags.append("unresolved_allele_alias")
    if observation.status == "ambiguous_reconstruction" and not (iupac or alias):
        flags.append("producer_status_unresolved")
    flags += [flag for key, flag in _COUNT_FLAGS if row.get(key)]
    return flags
