"""Allele length detection from samtools idxstats output.

Peak contig selection uses two alignment quality metrics extracted from the
BAM file:

- **Alignment score (AS):** Higher AS means the read aligns better to the
  contig.  Reads from a 60-repeat allele produce the highest AS when
  aligned to the contig whose length matches the allele (contig_51 for
  51 canonical X repeats + 9 fixed = 60 total).

- **Indel length:** Lower mean indel length indicates a better length
  match between read and reference.  Reads aligned to a contig that is
  too short or too long accumulate large insertions or deletions in the
  CIGAR string.

Both metrics independently identify the correct contig in testing.  We
use AS as the primary selector because it integrates all alignment
factors (matches, mismatches, gaps) into a single score.  Mean indel
length is reported alongside for transparency.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TypedDict

from muc_one_span.run_status import InsufficientEvidenceError
from muc_one_span.settings import DEFAULT_SETTINGS, AlleleSelectionSettings, ReferenceLayoutSettings
from muc_one_span.tools import run_tool_iter

# TypedDicts below document the expected structure of return values.
# Functions return plain dicts for mypy compatibility; these types are
# available for callers who want to annotate their own code.


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


logger = logging.getLogger(__name__)

# Number of fixed repeat units in the ladder reference (pre-repeats 1-5 + after-repeats 6-9).
# Each contig_N has N canonical X repeats plus these fixed repeats,
# so total allele length = N + PRE_AFTER_REPEAT_COUNT.
PRE_AFTER_REPEAT_COUNT = DEFAULT_SETTINGS.reference_layout.fixed_repeat_count


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


def _parse_cigar_indel_bp(cigar: str) -> int:
    """Sum of insertion and deletion bases from a CIGAR string."""
    total = 0
    for m in re.finditer(r"(\d+)([ID])", cigar):
        total += int(m.group(1))
    return total


def refine_peak_contig(
    bam_path: Path,
    cluster_contigs: list[str],
) -> dict:
    """Select the best contig from a cluster using alignment quality metrics.

    Scans all reads mapped to the cluster contigs and computes per-contig
    mean alignment score (AS tag) and mean indel length (from CIGAR).
    The contig with the **highest mean AS** is selected as the best match.

    Args:
        bam_path: Path to the ladder mapping BAM (indexed).
        cluster_contigs: List of contig names in the cluster
            (e.g. ``["contig_48", ..., "contig_54"]``).

    Returns:
        Dictionary with:

        - ``best_contig`` (str): name of the best-matching contig
        - ``metrics`` (dict): per-contig ``{mean_as, mean_indel_bp, reads}``
    """
    # Accumulate per-contig stats
    contig_stats: dict[str, dict] = {
        c: {
            "as_sum": 0,
            "indel_sum": 0,
            "count": 0,
            "primary": 0,
            "secondary": 0,
            "supplementary": 0,
        }
        for c in cluster_contigs
    }

    for line in run_tool_iter(["samtools", "view", str(bam_path), *cluster_contigs]):
        line = line.strip()
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue

        contig = fields[2]
        if contig not in contig_stats:
            continue

        flag = int(fields[1])
        if flag & 256:
            contig_stats[contig]["secondary"] += 1
        if flag & 2048:
            contig_stats[contig]["supplementary"] += 1
        if not flag & (4 | 256 | 2048):
            contig_stats[contig]["primary"] += 1
        cigar = fields[5]
        contig_stats[contig]["indel_sum"] += _parse_cigar_indel_bp(cigar)
        contig_stats[contig]["count"] += 1

        # Parse AS tag
        for tag in fields[11:]:
            if tag.startswith("AS:i:"):
                contig_stats[contig]["as_sum"] += int(tag[5:])
                break

    # Compute means and pick best
    metrics: dict[str, dict] = {}
    best_contig = cluster_contigs[0]
    best_as = -1.0

    for contig, stats in contig_stats.items():
        n = stats["count"]
        if n == 0:
            continue
        mean_as = stats["as_sum"] / n
        mean_indel = stats["indel_sum"] / n
        metrics[contig] = {
            "mean_as": round(mean_as, 1),
            "mean_indel_bp": round(mean_indel, 1),
            "reads": n,
            "primary_alignment_records": stats["primary"],
            "secondary_alignment_records": stats["secondary"],
            "supplementary_alignment_records": stats["supplementary"],
        }
        if mean_as > best_as:
            best_as = mean_as
            best_contig = contig

    logger.debug("Refined peak contig: %s (AS=%.1f)", best_contig, best_as)
    return {"best_contig": best_contig, "metrics": metrics}


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


def _split_cluster_by_indel(
    bam_path: Path,
    cluster: dict,
    *,
    settings: AlleleSelectionSettings | None = None,
) -> list[dict] | None:
    """Attempt to split a single cluster into two alleles using indel valleys.

    Reads from a short allele mapped to a long contig (or vice versa)
    accumulate large indels in CIGAR.  Reads mapped to the correct-length
    contig have near-zero indels.  By finding the two local minima in
    per-contig mean indel length, we can resolve close alleles that
    gap-based clustering merges into one cluster.

    Returns two sub-clusters if a clear split is found, or None when two strictly distinct valleys are not established.
    Failure to split does not establish sequence homozygosity.
    """
    settings = settings or DEFAULT_SETTINGS.allele_selection
    contig_names = [f"contig_{c}" for c, _ in cluster["contigs"]]

    # Compute per-contig mean indel bp
    contig_stats: dict[int, dict] = {}
    for c, _ in cluster["contigs"]:
        contig_stats[c] = {"indel_sum": 0, "count": 0}

    for line in run_tool_iter(["samtools", "view", str(bam_path), *contig_names]):
        line = line.strip()
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        contig_name = fields[2]
        m = re.search(r"_(\d+)$", contig_name)
        if not m:
            continue
        c = int(m.group(1))
        if c not in contig_stats:
            continue
        contig_stats[c]["indel_sum"] += _parse_cigar_indel_bp(fields[5])
        contig_stats[c]["count"] += 1

    # Build mean-indel series (only contigs with reads)
    indel_series: list[tuple[int, float]] = []
    for c in sorted(contig_stats):
        n = contig_stats[c]["count"]
        if n == 0:
            continue
        indel_series.append((c, contig_stats[c]["indel_sum"] / n))

    if len(indel_series) < settings.valley_min_points:
        return None

    # Find local minima (valleys) in mean indel
    valleys: list[tuple[int, float]] = []
    for i in range(len(indel_series)):
        c, val = indel_series[i]
        left = indel_series[i - 1][1] if i > 0 else float("inf")
        right = indel_series[i + 1][1] if i < len(indel_series) - 1 else float("inf")
        if val < left and val < right:
            valleys.append((c, val))

    logger.debug("Indel valley splitting: found %d valleys", len(valleys))
    if len(valleys) < 2:
        return None

    # Take the two deepest valleys (lowest mean indel)
    valleys.sort(key=lambda x: x[1])
    best_two = sorted(valleys[:2], key=lambda x: x[0])
    v1, v2 = best_two[0][0], best_two[1][0]

    # The split point is the midpoint between the two valleys
    split = (v1 + v2) // 2

    # Verify the valleys are meaningfully separated (at least 3 contigs apart)
    if abs(v2 - v1) < settings.valley_min_separation:
        return None

    # Split cluster contigs into two sub-clusters
    contigs_dict = dict(cluster["contigs"])
    sub1 = [(c, contigs_dict[c]) for c in sorted(contigs_dict) if c <= split]
    sub2 = [(c, contigs_dict[c]) for c in sorted(contigs_dict) if c > split]

    if not sub1 or not sub2:
        return None

    def _make_sub_cluster(contigs: list[tuple[int, int]]) -> dict:
        total = sum(r for _, r in contigs)
        center = sum(c * r for c, r in contigs) / total
        return {"center": round(center), "total_reads": total, "contigs": contigs}

    return [_make_sub_cluster(sub1), _make_sub_cluster(sub2)]


def _build_allele_info(
    cluster: dict,
    best_contig: str | None = None,
    *,
    settings: AlleleSelectionSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
) -> dict:
    """Build allele info dict from a cluster.

    Args:
        cluster: Cluster dict from _find_clusters.
        best_contig: Contig name selected by refine_peak_contig.
            If None, falls back to the weighted center.
    """
    settings = settings or DEFAULT_SETTINGS.allele_selection
    layout = reference_layout or DEFAULT_SETTINGS.reference_layout
    canonical = cluster["center"]
    contig_name = best_contig if best_contig is not None else f"contig_{canonical}"

    # When refine_peak_contig has identified a specific best contig,
    # derive canonical_repeats from its name rather than the cluster
    # center.  ONT reads produce wider distributions that can shift
    # the weighted center by ±1 vs the AS-refined best contig.
    # However, at longer allele lengths the AS metric itself can be
    # unreliable for ONT (off by 2+), so only trust the refinement
    # when it agrees within the configured shift (default ±1).
    if best_contig is not None:
        match = re.search(r"_(\d+)$", best_contig)
        if match:
            refined_canonical = int(match.group(1))
            if abs(refined_canonical - canonical) <= settings.refinement_max_shift:
                canonical = refined_canonical

    return {
        "length": canonical + layout.fixed_repeat_count,
        "fixed_repeat_count": layout.fixed_repeat_count,
        "reads": cluster["total_reads"],
        "alignment_records": cluster["total_reads"],
        "support_basis": "alignment_records",
        "molecule_count": None,
        "reference_length": int(contig_name.rsplit("_", 1)[1]) + layout.fixed_repeat_count,
        "canonical_repeats": canonical,
        "contig_name": contig_name,
        "cluster_contigs": [f"contig_{c}" for c, _ in cluster["contigs"]],
    }


def _length_selection_evidence(
    counts: dict[int, int],
    min_coverage: int,
    bam_path: Path | None,
    unselected_clusters: list[dict],
) -> dict:
    """Describe evidence omitted by the existing length-selection decisions."""
    excluded = sorted(
        (repeat, count) for repeat, count in counts.items() if 0 < count < min_coverage
    )
    primary_counts: dict[str, int] | None = None
    if excluded and bam_path is not None and bam_path.exists():
        contig_names = [f"contig_{repeat}" for repeat, _ in excluded]
        primary_counts = dict.fromkeys(contig_names, 0)
        for line in run_tool_iter(["samtools", "view", str(bam_path), *contig_names]):
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


def detect_alleles(
    counts: dict[int, int],
    min_coverage: int = DEFAULT_SETTINGS.run.min_coverage,
    bam_path: Path | None = None,
    *,
    settings: AlleleSelectionSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
) -> dict:
    """Detect allele lengths from read count distribution across ladder contigs.

    Finds two peak clusters in the read distribution. Each cluster represents
    one allele. Reports the total allele length (canonical repeats + configured
    fixed pre/after repeats) and the contig names needed for downstream processing.

    If *bam_path* is provided, the best contig within each cluster is refined
    using alignment scores (AS) and indel lengths from the BAM.  Without a
    BAM, the weighted center of the cluster is used as a fallback.

    Args:
        counts: Canonical repeat count -> mapped reads mapping from
            :func:`parse_idxstats`.
        min_coverage: Minimum mapped reads to include a contig.
        settings: Optional separation and refinement parameters.
        reference_layout: Selected fixed repeats used by the reference ladder.
        bam_path: Optional path to the indexed ladder mapping BAM.
            When provided, enables alignment-quality-based peak refinement.

    Returns:
        Dictionary with keys ``allele_1``, ``allele_2``, and ``homozygous``.
        Each allele has:

        - ``length`` (int): total repeat units including pre/after
        - ``reads`` (int): legacy count of alignment records across the cluster
        - ``alignment_records`` (int): explicit name for that same fit count
        - ``molecule_count``: unavailable without validated molecular identity
        - ``reference_length`` (int): total repeats in the selected reference
        - ``fixed_repeat_count`` (int): selected fixed pre/after repeat count
        - ``canonical_repeats`` (int): number of canonical X repeats
        - ``contig_name`` (str): best-matching contig (e.g. ``"contig_51"``)
        - ``cluster_contigs`` (list[str]): all contig names in the cluster
        - ``length_selection_evidence``: diagnostic accounting for sub-threshold
          contigs and passing clusters not selected among the first two

    Raises:
        ValueError: If no contig meets the minimum coverage threshold.
    """
    settings = settings or DEFAULT_SETTINGS.allele_selection
    # Handle mixed-key dicts (legacy compat): only use integer keys
    int_counts = {k: v for k, v in counts.items() if isinstance(k, int)}

    clusters = _find_clusters(int_counts, min_coverage, settings.min_gap)
    logger.info("Detected %d peak(s) from %d contigs", len(clusters), len(int_counts))

    if not clusters:
        max_observed = max(int_counts.values()) if int_counts else 0
        raise InsufficientEvidenceError(
            f"No contig has >= {min_coverage} mapped reads (minimum coverage). "
            f"Max observed: {max_observed} reads."
        )

    fit_metrics: dict[str, dict] = {}

    # Refine peak contig selection using alignment quality if BAM available
    def _get_best_contig(cluster: dict) -> str | None:
        if bam_path is None:
            return None
        contig_names = [f"contig_{c}" for c, _ in cluster["contigs"]]
        refined = refine_peak_contig(bam_path, contig_names)
        fit_metrics.update(refined["metrics"])
        best: str = refined["best_contig"]
        return best

    # If only one cluster found but BAM is available, try indel-valley splitting
    if len(clusters) == 1 and bam_path is not None:
        sub_clusters = _split_cluster_by_indel(bam_path, clusters[0], settings=settings)
        if sub_clusters is not None:
            clusters = sub_clusters
            clusters.sort(key=lambda x: x["total_reads"], reverse=True)

    selection_evidence = _length_selection_evidence(
        int_counts, min_coverage, bam_path, clusters[2:]
    )

    def _with_support(cluster: dict) -> dict:
        info = _build_allele_info(
            cluster, _get_best_contig(cluster), settings=settings, reference_layout=reference_layout
        )
        info["primary_alignment_records"] = (
            sum(
                fit_metrics.get(c, {}).get("primary_alignment_records", 0)
                for c in info["cluster_contigs"]
            )
            if bam_path is not None
            else None
        )
        info["fit_metrics"] = {
            c: fit_metrics[c] for c in info["cluster_contigs"] if c in fit_metrics
        }
        info["length_selection_evidence"] = selection_evidence
        return info

    allele_1 = _with_support(clusters[0])

    if len(clusters) < 2:
        return {
            "allele_1": allele_1,
            "allele_2": {**allele_1, "candidate_duplicate_of": "allele_1"},
            "observed_length_candidates": 1,
            "allele_multiplicity_status": "unresolved",
            "homozygous": False,
            "same_length": True,
        }

    allele_2 = _with_support(clusters[1])

    if allele_1["length"] == allele_2["length"]:
        allele_1["reads"] += allele_2["reads"]
        allele_1["alignment_records"] = allele_1["reads"]
        allele_1["cluster_contigs"] = sorted(
            set(allele_1["cluster_contigs"] + allele_2["cluster_contigs"])
        )
        allele_1["fit_metrics"].update(allele_2["fit_metrics"])
        if bam_path is not None:
            allele_1["primary_alignment_records"] = sum(
                m.get("primary_alignment_records", 0) for m in allele_1["fit_metrics"].values()
            )
        return {
            "allele_1": allele_1,
            "allele_2": {**allele_1, "candidate_duplicate_of": "allele_1"},
            "observed_length_candidates": 1,
            "allele_multiplicity_status": "unresolved",
            "homozygous": False,
            "same_length": True,
        }

    return {
        "allele_1": allele_1,
        "allele_2": allele_2,
        "homozygous": False,
        "same_length": False,
        "observed_length_candidates": len(clusters),
        "allele_multiplicity_status": "two_selected_candidates",
    }
