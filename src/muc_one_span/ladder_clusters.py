"""Pure clustering helpers for allele detection from samtools idxstats output."""

from __future__ import annotations

import re
from typing import TypedDict

from muc_one_span.settings import DEFAULT_SETTINGS


class AlleleInfo(TypedDict):
    """Information about a single detected allele."""

    length: int
    reads: int
    canonical_repeats: int
    contig_name: str
    cluster_contigs: list[str]


class AlleleResult(TypedDict):
    """Result of allele detection for a sample."""

    allele_1: AlleleInfo
    allele_2: AlleleInfo
    homozygous: bool
    same_length: bool


def parse_idxstats(idxstats_output: str) -> dict[int, int]:
    """Parse samtools idxstats output into repeat_count -> read_count mapping.

    Expects contig names like 'contig_60' where 60 is the number of
    canonical X repeats in that contig.

    Args:
        idxstats_output: Raw text output from ``samtools idxstats``.

    Returns:
        Dictionary mapping canonical repeat count (int) to mapped alignment
        record count (int), including secondary alignments.  The '*' unmapped line is excluded.
    """
    counts: dict[int, int] = {}

    for line in idxstats_output.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue

        contig_name = parts[0]
        if contig_name == "*":
            continue

        mapped_reads = int(parts[2])

        match = re.search(r"_(\d+)$", contig_name)
        if match:
            repeat_count = int(match.group(1))
            counts[repeat_count] = mapped_reads

    return counts


def _find_clusters(
    counts: dict[int, int],
    min_coverage: int,
    min_gap: int = DEFAULT_SETTINGS.allele_selection.min_gap,
) -> list[dict]:
    """Identify read-count clusters in the contig distribution.

    Groups contigs that are within ``min_gap`` of each other into clusters,
    then computes the weighted center and total reads for each.

    Args:
        counts: Canonical repeat count -> mapped reads mapping.
        min_coverage: Minimum reads for a contig to be included.
        min_gap: Minimum gap between contigs to start a new cluster.

    Returns:
        List of cluster dicts sorted by total_reads descending.
        Each dict has keys: center (int), total_reads (int),
        contigs (list of (repeat_count, reads) tuples).
    """
    passing = sorted(
        [(k, v) for k, v in counts.items() if v >= min_coverage],
        key=lambda x: x[0],
    )

    if not passing:
        return []

    # Group into clusters by proximity
    clusters: list[list[tuple[int, int]]] = []
    current_cluster: list[tuple[int, int]] = [passing[0]]

    for i in range(1, len(passing)):
        if passing[i][0] - passing[i - 1][0] >= min_gap:
            clusters.append(current_cluster)
            current_cluster = [passing[i]]
        else:
            current_cluster.append(passing[i])
    clusters.append(current_cluster)

    # Compute weighted center and total reads for each cluster
    result: list[dict] = []
    for cluster in clusters:
        total_reads = sum(reads for _, reads in cluster)
        weighted_center = sum(pos * reads for pos, reads in cluster) / total_reads
        result.append(
            {
                "center": round(weighted_center),
                "total_reads": total_reads,
                "contigs": cluster,
            }
        )

    result.sort(key=lambda x: x["total_reads"], reverse=True)
    return result
