"""Typed result schemas for repeat and sequence classification."""

from __future__ import annotations

from typing import TypedDict

# TypedDicts below document the expected structure of return values.
# Functions return plain dicts for mypy compatibility; these types are
# available for callers who want to annotate their own code.


class RepeatDifference(TypedDict):
    """A single difference between a repeat sequence and its closest reference."""

    pos: int
    ref: str
    alt: str
    type: str


class RepeatClassification(TypedDict, total=False):
    """Classification result for a single repeat unit."""

    type: str
    match: str
    confidence: float
    closest_match: str
    edit_distance: int | float
    identity_pct: float
    differences: list[RepeatDifference]
    classification: str
    frameshift: bool
    mutation_name: str
    parent_repeat: str
    index: int


class MutationDetected(TypedDict, total=False):
    """A mutation detected during classification."""

    repeat_index: int
    closest_type: str
    mutation_name: str
    template_match: bool
    frameshift: bool
    differences: list[RepeatDifference]
    vcf_support: bool
    vcf_qual: float
    boundary: bool


class SequenceClassification(TypedDict):
    """Classification result for a full consensus sequence."""

    structure: str
    repeats: list[RepeatClassification]
    mutations_detected: list[MutationDetected]
    cumulative_offset: int
    allele_confidence: float
    exact_match_pct: float
