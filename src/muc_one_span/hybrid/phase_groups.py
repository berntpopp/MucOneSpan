"""Task 15k: exclude a linked-site group explained by a low-accuracy read subset.

In a peak of a two-peak model each allele already has its own length peak. On
simulated HiFi amplicons, a subset of lower-quality reads of one allele can share two
systematic errors a few bases apart (a G read as C twice inside one unit, say); the
errors are linked through the same reads, so the peak is split into two groups and the
sample carries a third allele group (``unresolved_max_alleles``). Task 15j keeps that
split, because merging the poor reads back would put them into the allele consensus.

This rule keeps the split too, and never lets the poor group join an allele consensus.
It only stops counting that group as a further allele when every condition holds:

1. it is the smaller group of the split;
2. its reads have lower mean base quality than the other group's (one-sided
   rank-sum test at ``phase_quality_alpha``, as in ``phase_quality``);
3. among the ``phase_quality_keep_frac`` of the split's reads with the highest mean
   base quality, its share is significantly below ``het_af_min`` (exact binomial
   upper bound at ``phase_quality_af_alpha``): it is not a real minority haplotype
   carried by good reads;
4. it carries no distinct length: its median read length is within
   ``peak_min_separation_units`` repeat units of the other group's;
5. it carries no event: its draft differs from the other group's draft by
   substitutions only (same length, and the edit distance equals the number of
   mismatched positions), so it holds no insertion or deletion and hence no
   frameshift-class event.

It fails closed like ``phase_quality``: too few high-quality reads, a read without
base qualities, or no quality variance keep the group (and the unresolved status).
The excluded reads are counted as spanning reads assigned to no allele, so the
``max_unassigned_spanning_fraction`` guard still sees them.
"""

from __future__ import annotations

import math
import statistics
from typing import Any

from muc_one_span.hybrid.align import edit_distance
from muc_one_span.hybrid.allele_fields import PLOIDY
from muc_one_span.hybrid.phase_quality import (
    NO_QUALITY,
    P_DIGITS,
    binomial_lower_tail,
    quality_floor_reads,
    rank_sum_p_lower,
)
from muc_one_span.hybrid.phase_sites import AF_DECIMALS
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings


def substitutions_only(a: str, b: str) -> bool:
    """True when ``a`` and ``b`` have the same length and differ only by mismatches."""
    if len(a) != len(b):
        return False
    return edit_distance(a, b) == sum(x != y for x, y in zip(a, b, strict=True))


def explained_group(
    groups: list[list[SpanRead]], drafts: list[str], s: HybridSettings, unit_bp: int
) -> tuple[int, dict[str, Any]] | None:
    """(index, evidence) of a two-group split's group explained by poor reads, or None.

    ``groups`` and ``drafts`` are the split's two groups and their POA drafts, in the
    same order; ``unit_bp`` is the repeat-unit length (from the dictionary).
    """
    if not s.phase_quality_group_exclusion or not len(groups) == len(drafts) == PLOIDY:
        return None
    small, large = sorted(range(len(groups)), key=lambda i: len(groups[i]))
    quals = [[m.mean_q for m in g] for g in groups]
    if any(q <= NO_QUALITY for g in quals for q in g):
        return None
    p = rank_sum_p_lower(quals[small], quals[large])
    if p >= s.phase_quality_alpha:
        return None
    ranked = sorted(
        [(q, True) for q in quals[small]] + [(q, False) for q in quals[large]],
        key=lambda item: -item[0],
    )
    kept = ranked[: math.ceil(s.phase_quality_keep_frac * len(ranked))]
    if len(kept) < quality_floor_reads(s):
        return None
    k = sum(is_small for _q, is_small in kept)
    tail = binomial_lower_tail(k, len(kept), s.het_af_min)
    if tail >= s.phase_quality_af_alpha:
        return None
    lengths = [statistics.median(m.length for m in g) for g in groups]
    if abs(lengths[small] - lengths[large]) >= s.peak_min_separation_units * unit_bp:
        return None
    if not substitutions_only(drafts[small], drafts[large]):
        return None
    return small, {
        "spanning_reads": len(groups[small]),
        "median_units": round(lengths[small] / unit_bp),
        "share_high_quality": round(k / len(kept), AF_DECIMALS),
        "quality_p": float(f"{p:.{P_DIGITS}g}"),
        "af_bound_p": float(f"{tail:.{P_DIGITS}g}"),
        "draft_substitutions": sum(
            x != y for x, y in zip(drafts[small], drafts[large], strict=True)
        ),
    }
