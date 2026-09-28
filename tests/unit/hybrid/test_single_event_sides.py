"""The single-event share gate bounds both sides of the split.

The gate's minor allele is the site's non-draft allele. When the POA draft follows an
artefact (a homozygous normal with a delinsAT-shaped artefact in 40% of its reads was
drafted with the artefact's C4 run), the "minor" allele is the wild type at 60%: its
bound cleared ``phase_single_event_min_share`` and the peak was split, and the 40%
artefact group was called PATHOGENIC. Both the minor share and the major share must
reach the floor, so the smaller group of any split is at least the floor.
"""

from __future__ import annotations

from typing import Any

from muc_one_span.hybrid import single_event
from tests.unit.hybrid import test_single_event as base

S = base.S
SITE: dict[str, Any] = {"site": ("col", 5), "major": "C", "minor": "A"}


def _feats(n_minor: int, n_major: int) -> tuple[list[dict[Any, Any]], list[str]]:
    feats = [{SITE["site"]: SITE["minor"]}] * n_minor + [{SITE["site"]: SITE["major"]}] * n_major
    strands = ["+", "-"] * (len(feats) // 2)
    return feats, strands


def test_both_shares_are_bounded() -> None:
    total = S.phase_single_event_bound_reads
    minor = round((1 - S.phase_single_event_min_share) * total)
    feats, strands = _feats(minor, total - minor)
    minor_bound, major_bound = single_event._share_bounds(feats, strands, SITE, None, S)
    # The draft follows the smaller allele: its "minor" side clears the floor ...
    assert minor_bound >= S.phase_single_event_min_share
    # ... but the draft's own allele, at exactly the floor, cannot be bounded above it.
    assert major_bound < S.phase_single_event_min_share


def test_a_balanced_site_passes_on_both_sides() -> None:
    half = S.phase_single_event_bound_reads // 2
    feats, strands = _feats(half, half)
    bounds = single_event._share_bounds(feats, strands, SITE, None, S)
    assert min(bounds) >= S.phase_single_event_min_share
