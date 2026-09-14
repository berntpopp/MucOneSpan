"""Offline evaluation records; positions/counts include every VNTR repeat unit."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from muc_one_span.settings import DEFAULT_SETTINGS


@dataclass(frozen=True)
class Event:
    """Exact annotation identity; this is not normalized nucleotide equivalence."""

    repeat_index: int | None
    parent: str | None
    name: str | None
    frameshift: bool = field(default=False, compare=False)
    template_match: bool = field(default=False, compare=False)
    vcf_support: bool = field(default=False, compare=False)
    support_status: str = field(default="unknown", compare=False)

    @property
    def supported(self) -> bool:
        """Return whether exact template and sequence-projection evidence agree."""
        return self.legacy_supported and self.support_status == "exact_sequence_concordance"

    @property
    def legacy_supported(self) -> bool:
        """Return the historical boolean support rule without projection proof."""
        return self.frameshift and self.template_match and self.vcf_support


@dataclass(frozen=True)
class TruthHaplotype:
    name: str
    sequence: str
    structure: tuple[str, ...]
    events: tuple[Event, ...] = ()


@dataclass(frozen=True)
class TruthSample:
    name: str
    haplotypes: tuple[TruthHaplotype, ...]
    provenance: dict[str, Any] = field(default_factory=dict)
    hashes: dict[str, str] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class PredictedAllele:
    name: str
    sequence: str
    structure: tuple[str, ...]
    length: int
    canonical_repeats: int
    events: tuple[Event, ...] = ()
    reconstruction_status: str = "unknown"
    sequence_source: str | None = None
    independent_haplotype_evidence: bool = False
    phase_status: str = "unknown"
    genotype_status: str = "unknown"
    sequence_identity_status: str = "unknown"
    fixed_repeat_count: int = DEFAULT_SETTINGS.reference_layout.fixed_repeat_count


@dataclass(frozen=True)
class RunObservation:
    status: str
    predictions: tuple[PredictedAllele, ...] = ()
    run_record: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    error: str | None = None


@dataclass(frozen=True)
class Assignment:
    cost: int
    alternatives: tuple[tuple[tuple[int, int], ...], ...]
    distance_matrix: tuple[tuple[int, ...], ...]
