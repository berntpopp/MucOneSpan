# src/muc_one_span/read_dominance.py
"""Read-dominance rescoring and candidate validation for allele length inference."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from muc_one_span.settings import DEFAULT_SETTINGS
from muc_one_span.tools import run_tool_iter

logger = logging.getLogger(__name__)

# Platform defaults, sourced from ``allele_selection`` (single source of truth).
# Score margins derive from minimap2 two-piece affine gap minimum g: delta = 0.5 * g.
# map-hifi: -O6,26 -E2,1 -> min(6+120, 26+60) = 86 -> delta = 43
# lr:hq:    -O4,24 -E2,1 -> min(4+120, 24+60) = 84 -> delta = 42
_SELECTION = DEFAULT_SETTINGS.allele_selection
DEFAULT_DELTA_AS = {
    "hifi": _SELECTION.score_margin_hifi,
    "ont": _SELECTION.score_margin_ont,
    "default": _SELECTION.score_margin_hifi,
}
DEFAULT_MIN_DOMINANT_READS = {
    "hifi": _SELECTION.min_dominant_reads_hifi,
    "ont": _SELECTION.min_dominant_reads_ont,
    "default": _SELECTION.min_dominant_reads_hifi,
}
DEFAULT_MIN_RATIO = _SELECTION.min_dominance_ratio


@dataclass(frozen=True)
class DominanceScore:
    """Pairwise dominance scores between two candidate contigs."""

    c1: str
    c2: str
    d1_count: int  # Reads strictly preferring c1 by > delta
    d2_count: int  # Reads strictly preferring c2 by > delta
    ambiguous_count: int  # Reads within delta margin
    total_scored_reads: int
    is_valid_second_allele: bool
    rejection_reason: str | None = None


def extract_read_scores_for_contigs(
    bam_path: Path,
    contig_names: list[str],
    *,
    run_tool_iter_func: Any = None,
) -> dict[str, dict[str, int]]:
    """Extract alignment score (AS) for each read across specified contigs.

    Args:
        bam_path: Path to coordinate-sorted, indexed BAM file.
        contig_names: List of contig names to extract (e.g. ['contig_51', 'contig_54']).
        run_tool_iter_func: Optional runner function (defaults to run_tool_iter).

    Returns:
        Mapping of read_qname -> {contig_name: alignment_score}.
    """
    scores: dict[str, dict[str, int]] = {}
    valid_contigs = set(contig_names)
    if not valid_contigs or not bam_path.exists():
        return scores

    runner = run_tool_iter_func or run_tool_iter
    cmd = ["samtools", "view", "-F", "0x800", str(bam_path), *contig_names]
    for line in runner(cmd):
        line = line.strip()
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        qname = fields[0]
        rname = fields[2]
        if rname not in valid_contigs:
            continue

        as_score: int | None = None
        for tag in fields[11:]:
            if tag.startswith("AS:i:"):
                try:
                    as_score = int(tag[5:])
                except ValueError:
                    as_score = None
                break

        if as_score is not None:
            if qname not in scores:
                scores[qname] = {}
            # Retain highest score if multiple alignments to same contig
            if rname not in scores[qname] or as_score > scores[qname][rname]:
                scores[qname][rname] = as_score

    return scores


def evaluate_candidate_pair_dominance(
    read_scores: dict[str, dict[str, int]],
    c1: str,
    c2: str,
    platform: str = "hifi",
    *,
    delta: int | None = None,
    min_dominant_reads: int | None = None,
    min_ratio: float = DEFAULT_MIN_RATIO,
    c2_primary_records: int | None = None,
    close_candidate_repeats: int | None = None,
    zero_primary_extra_reads: int | None = None,
) -> DominanceScore:
    """Evaluate whether candidate c2 has genuine read-dominance over c1.

    Args:
        read_scores: Mapping from extract_read_scores_for_contigs.
        c1: Dominant or higher-coverage candidate contig name.
        c2: Secondary or alternative candidate contig name.
        platform: 'hifi' or 'ont'.
        delta: Score difference threshold; None uses platform default.
        min_dominant_reads: Minimum count of reads strictly preferring c2; None uses platform default.
        min_ratio: Minimum ratio |D2| / (|D1| + |D2|).
        c2_primary_records: Optional count of primary alignments on c2.
        close_candidate_repeats: Repeat difference below which a read aligned to only
            one candidate is ambiguous; None uses
            ``allele_selection.dominance_close_candidate_repeats``.
        zero_primary_extra_reads: Extra dominant reads required when c2 has no primary
            alignment; None uses ``allele_selection.dominance_zero_primary_extra_reads``.

    Returns:
        DominanceScore with decision and detailed counts.
    """
    effective_delta = (
        delta if delta is not None else DEFAULT_DELTA_AS.get(platform, DEFAULT_DELTA_AS["default"])
    )
    effective_min_reads = (
        min_dominant_reads
        if min_dominant_reads is not None
        else DEFAULT_MIN_DOMINANT_READS.get(platform, DEFAULT_MIN_DOMINANT_READS["default"])
    )

    close_repeats = (
        close_candidate_repeats
        if close_candidate_repeats is not None
        else _SELECTION.dominance_close_candidate_repeats
    )
    extra_reads = (
        zero_primary_extra_reads
        if zero_primary_extra_reads is not None
        else _SELECTION.dominance_zero_primary_extra_reads
    )

    d1 = 0
    d2 = 0
    ambiguous = 0

    try:
        n1 = int(c1.split("_")[-1])
        n2 = int(c2.split("_")[-1])
        delta_repeats = abs(n1 - n2)
    except (ValueError, IndexError):
        delta_repeats = 0

    for _qname, c_scores in read_scores.items():
        s1 = c_scores.get(c1)
        s2 = c_scores.get(c2)
        if s1 is None and s2 is None:
            continue
        elif s1 is not None and s2 is None:
            if delta_repeats >= close_repeats:
                d1 += 1
            else:
                ambiguous += 1
        elif s2 is not None and s1 is None:
            if delta_repeats >= close_repeats:
                d2 += 1
            else:
                # Close candidate: secondary suppression or stutter cannot be excluded
                ambiguous += 1
        else:
            assert s1 is not None and s2 is not None
            diff = s2 - s1
            if diff > effective_delta:
                d2 += 1
            elif diff < -effective_delta:
                d1 += 1
            else:
                ambiguous += 1

    total_scored = d1 + d2 + ambiguous

    # Rejection conditions
    if d2 < effective_min_reads:
        return DominanceScore(
            c1=c1,
            c2=c2,
            d1_count=d1,
            d2_count=d2,
            ambiguous_count=ambiguous,
            total_scored_reads=total_scored,
            is_valid_second_allele=False,
            rejection_reason=f"Insufficient dominant reads: {d2} < {effective_min_reads}",
        )

    decisive_total = d1 + d2
    if decisive_total > 0 and (d2 / decisive_total) < min_ratio:
        return DominanceScore(
            c1=c1,
            c2=c2,
            d1_count=d1,
            d2_count=d2,
            ambiguous_count=ambiguous,
            total_scored_reads=total_scored,
            is_valid_second_allele=False,
            rejection_reason=f"Dominance ratio below threshold: {d2}/{decisive_total} < {min_ratio}",
        )

    if (
        c2_primary_records is not None
        and c2_primary_records == 0
        and d2 < (effective_min_reads + extra_reads)
    ):
        return DominanceScore(
            c1=c1,
            c2=c2,
            d1_count=d1,
            d2_count=d2,
            ambiguous_count=ambiguous,
            total_scored_reads=total_scored,
            is_valid_second_allele=False,
            rejection_reason="No primary alignments on secondary candidate and marginal dominance",
        )

    return DominanceScore(
        c1=c1,
        c2=c2,
        d1_count=d1,
        d2_count=d2,
        ambiguous_count=ambiguous,
        total_scored_reads=total_scored,
        is_valid_second_allele=True,
    )
