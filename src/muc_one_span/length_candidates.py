# src/muc_one_span/length_candidates.py
"""Candidate length discovery and refinement for MUC1 VNTR alleles."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.read_dominance import (
    DominanceScore,
    evaluate_candidate_pair_dominance,
    extract_read_scores_for_contigs,
)
from muc_one_span.settings import (
    DEFAULT_SETTINGS,
    AlleleSelectionSettings,
    ReferenceLayoutSettings,
)
from muc_one_span.tools import run_tool_iter

logger = logging.getLogger(__name__)

# Structural, not tunable: a read-length peak must dominate the bins up to this many
# bin widths away on each side (the former ``b +/- 5`` and ``b +/- 10`` neighbours).
PEAK_NEIGHBOUR_BINS = 2


@dataclass
class CandidateCluster:
    """Discovered candidate length cluster."""

    canonical_repeats: int
    contig_name: str
    total_reads: int
    primary_reads: int
    cluster_contigs: list[str]


def discover_length_candidates(
    counts: dict[int, int],
    primary_counts: dict[int, int] | None = None,
    min_coverage: int = 10,
    discovery_floor: int = 3,
) -> list[CandidateCluster]:
    """Discover candidate allele lengths from read count distribution.

    Identifies primary modes and potential minority candidates.

    Args:
        counts: Mapping of canonical repeat count -> total mapped reads.
        primary_counts: Optional mapping of canonical repeat count -> primary alignments.
        min_coverage: Threshold for standard coverage alleles.
        discovery_floor: Lower threshold for minority/asymmetric allele candidates.

    Returns:
        List of CandidateCluster objects sorted by evidence strength.
    """
    prim = primary_counts or {}
    candidates: list[CandidateCluster] = []

    # Find peaks with >= min_coverage
    eligible = {k: v for k, v in counts.items() if isinstance(k, int) and v >= discovery_floor}
    if not eligible:
        return []

    # Sort contigs by reads
    sorted_contigs = sorted(eligible.items(), key=lambda item: item[1], reverse=True)
    best_c, best_reads = sorted_contigs[0]

    primary_candidate = CandidateCluster(
        canonical_repeats=best_c,
        contig_name=f"contig_{best_c}",
        total_reads=best_reads,
        primary_reads=prim.get(best_c, 0),
        cluster_contigs=[f"contig_{best_c}"],
    )
    candidates.append(primary_candidate)

    # Search for a second distinct candidate
    # Candidate must be separated from best_c or have distinct primary support
    secondary_candidates: list[CandidateCluster] = []
    for c, r in sorted_contigs[1:]:
        if c == best_c:
            continue
        p = prim.get(c, 0)
        # To qualify as a secondary candidate:
        # Either has primary alignments >= 1, or has >= min_coverage, or is distant (|c - best_c| >= 5) with >= discovery_floor
        is_separated = abs(c - best_c) >= 5 and r >= discovery_floor
        has_primary = p >= 1 and r >= discovery_floor
        has_coverage = r >= min_coverage

        if is_separated or has_primary or has_coverage:
            secondary_candidates.append(
                CandidateCluster(
                    canonical_repeats=c,
                    contig_name=f"contig_{c}",
                    total_reads=r,
                    primary_reads=p,
                    cluster_contigs=[f"contig_{c}"],
                )
            )

    if secondary_candidates:
        # Sort secondary candidates by primary alignments then total reads
        secondary_candidates.sort(key=lambda x: (x.primary_reads, x.total_reads), reverse=True)
        candidates.append(secondary_candidates[0])

    return candidates


def validate_candidates_with_dominance(
    bam_path: Path | None,
    candidate_1: CandidateCluster,
    candidate_2: CandidateCluster | None,
    platform: str = "hifi",
    settings: AlleleSelectionSettings | None = None,
) -> tuple[bool, DominanceScore | None]:
    """Validate whether candidate_2 is supported by read-level dominance over candidate_1.

    Args:
        bam_path: Path to BAM file for reading alignment scores.
        candidate_1: Primary candidate cluster.
        candidate_2: Secondary candidate cluster, if proposed.
        platform: 'hifi' or 'ont'.
        settings: Allele selection settings.

    Returns:
        Tuple of (is_valid_second_allele, dominance_score_object).
    """
    if candidate_2 is None or candidate_1.canonical_repeats == candidate_2.canonical_repeats:
        return False, None

    if bam_path is None or not bam_path.exists():
        # Without BAM, cannot perform dominance rescoring; rely on clustering
        return True, None

    c1_name = candidate_1.contig_name
    c2_name = candidate_2.contig_name

    # Extract alignment scores for both contigs
    read_scores = extract_read_scores_for_contigs(bam_path, [c1_name, c2_name])
    if not read_scores:
        return False, None

    dom_score = evaluate_candidate_pair_dominance(
        read_scores,
        c1_name,
        c2_name,
        platform=platform,
        c2_primary_records=candidate_2.primary_reads,
    )

    return dom_score.is_valid_second_allele, dom_score


def _length_selection_evidence(
    counts: dict[int, int],
    min_coverage: int,
    bam_path: Path | None,
    unselected_clusters: list[dict],
    run_tool_iter_func: Any = None,
) -> dict:
    """Describe evidence omitted by the existing length-selection decisions."""
    runner = run_tool_iter_func or run_tool_iter
    excluded = sorted(
        (repeat, count) for repeat, count in counts.items() if 0 < count < min_coverage
    )
    primary_counts: dict[str, int] | None = None
    if excluded and bam_path is not None and bam_path.exists():
        contig_names = [f"contig_{repeat}" for repeat, _ in excluded]
        primary_counts = dict.fromkeys(contig_names, 0)
        for line in runner(["samtools", "view", str(bam_path), *contig_names]):
            fields = line.strip().split("\t")
            if len(fields) < 3 or fields[2] not in primary_counts:
                continue
            flag = int(fields[1])
            if not flag & (4 | 256 | 2048):
                primary_counts[fields[2]] += 1

    excluded_rows = [
        {
            "contig_name": f"contig_{repeat}",
            "alignment_records": count,
            "primary_alignment_records": (
                primary_counts[f"contig_{repeat}"] if primary_counts is not None else None
            ),
            "molecule_count": None,
        }
        for repeat, count in excluded
    ]
    if primary_counts is not None:
        excluded_primary: int | None = sum(primary_counts.values())
    else:
        excluded_primary = 0 if not excluded else None

    return {
        "minimum_coverage": min_coverage,
        "support_unit": "alignment_records_not_molecules",
        "molecule_count": None,
        "excluded_subthreshold_contigs": excluded_rows,
        "excluded_subthreshold_alignment_records": sum(count for _, count in excluded),
        "excluded_subthreshold_primary_alignment_records": excluded_primary,
        "unselected_passing_clusters": [
            {
                "center": cluster["center"],
                "alignment_records": cluster["total_reads"],
                "contigs": [
                    {"contig_name": f"contig_{repeat}", "alignment_records": count}
                    for repeat, count in cluster["contigs"]
                ],
            }
            for cluster in unselected_clusters
        ],
    }


def split_cluster_by_read_length(
    bam_path: Path,
    cluster: dict,
    *,
    platform: str = "hifi",
    min_reads: int | None = None,
    min_fraction: float | None = None,
    run_tool_iter_func: Any = None,
    settings: AlleleSelectionSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
    repeat_length_bp: int | None = None,
) -> list[dict] | None:
    """Split a single cluster into two alleles if physical read lengths are bimodal.

    In spanning reads (particularly PacBio HiFi amplicons), read length directly
    reflects VNTR repeat count (one repeat unit per ``repeat_length_bp``). Near-equal
    alleles (e.g. 40/41) merge into one cluster during gap-based contig clustering,
    but exhibit two distinct modes in physical read lengths.

    Thresholds come from ``settings`` (``allele_selection.read_length_split_*``);
    *min_reads* and *min_fraction* override their settings when given. The unit length
    defaults to the bundled repeat dictionary's ``repeat_length_bp`` and the fixed
    repeat count to *reference_layout*. A peak's bin count must be at least that of
    every bin up to :data:`PEAK_NEIGHBOUR_BINS` bin widths away on each side (a
    structural constant, not a setting).

    Returns two sub-clusters if two valid modes are detected, or None.

    Raises:
        ValueError: If twice the unit tolerance is not below the repeat unit length.
    """
    s = settings or DEFAULT_SETTINGS.allele_selection
    layout = reference_layout or DEFAULT_SETTINGS.reference_layout
    unit = (
        load_repeat_dictionary().repeat_length_bp if repeat_length_bp is None else repeat_length_bp
    )
    tolerance = s.read_length_split_unit_tolerance_bp
    if 2 * tolerance >= unit:
        raise ValueError(
            f"allele_selection.read_length_split_unit_tolerance_bp={tolerance} must be "
            f"below half the repeat unit length ({unit} bp)"
        )
    min_reads = s.read_length_split_min_reads if min_reads is None else min_reads
    min_fraction = s.read_length_split_min_fraction if min_fraction is None else min_fraction
    bin_bp = s.read_length_split_bin_bp

    runner = run_tool_iter_func or run_tool_iter
    contig_names = [f"contig_{c}" for c, _ in cluster["contigs"]]
    if not contig_names or not bam_path.exists():
        return None

    read_lengths: list[int] = []
    for line in runner(["samtools", "view", "-F", "2308", str(bam_path), *contig_names]):
        fields = line.strip().split("\t")
        if len(fields) >= 10 and fields[9] != "*":
            read_lengths.append(len(fields[9]))

    total_reads = len(read_lengths)
    if total_reads < min_reads * 2:
        return None

    bins: dict[int, int] = {}
    for length in read_lengths:
        b = (length // bin_bp) * bin_bp
        bins[b] = bins.get(b, 0) + 1

    peaks: list[tuple[int, int]] = []
    threshold = max(min_reads, int(total_reads * min_fraction))
    offsets = [k * bin_bp for k in range(1, PEAK_NEIGHBOUR_BINS + 1)]
    for b in sorted(bins):
        cnt = bins[b]
        if cnt >= threshold:
            left = max(bins.get(b - off, 0) for off in offsets)
            right = max(bins.get(b + off, 0) for off in offsets)
            if cnt >= left and cnt >= right:
                peaks.append((b, cnt))

    if len(peaks) < 2:
        return None

    peaks.sort(key=lambda x: x[1], reverse=True)
    p1, p2 = sorted([peaks[0][0], peaks[1][0]])
    delta = p2 - p1

    if not (s.read_length_split_min_delta_bp <= delta <= s.read_length_split_max_delta_bp):
        return None
    rem = delta % unit
    if not (rem <= tolerance or rem >= unit - tolerance):
        return None

    offset = s.read_length_split_offset_bp
    c1 = round((p1 - offset) / unit) - layout.fixed_repeat_count
    c2 = round((p2 - offset) / unit) - layout.fixed_repeat_count
    if c1 >= c2:
        return None

    contigs_dict = dict(cluster["contigs"])
    sub1 = [(c, contigs_dict[c]) for c in sorted(contigs_dict) if abs(c - c1) < abs(c - c2)]
    sub2 = [(c, contigs_dict[c]) for c in sorted(contigs_dict) if abs(c - c2) < abs(c - c1)]

    if not any(c == c1 for c, _ in sub1):
        sub1.append((c1, max(1, peaks[0][1])))
        sub1.sort(key=lambda x: x[0])
    if not any(c == c2 for c, _ in sub2):
        sub2.append((c2, max(1, peaks[1][1])))
        sub2.sort(key=lambda x: x[0])

    def _make_sub(contigs: list[tuple[int, int]], peak_c: int) -> dict:
        t = sum(r for _, r in contigs)
        return {
            "center": peak_c,
            "total_reads": t,
            "contigs": contigs,
            "split_diagnostics": {
                "splitter": "read_length",
                "peaks": [p1, p2],
                "delta": delta,
                "target_c": peak_c,
                "read_count": total_reads,
            },
        }

    return [_make_sub(sub1, c1), _make_sub(sub2, c2)]
