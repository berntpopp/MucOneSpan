"""Allele length detection from samtools idxstats output and alignment metrics."""

from __future__ import annotations

import logging
import re
from pathlib import Path

from muc_one_span.ladder_clusters import AlleleInfo as AlleleInfo
from muc_one_span.ladder_clusters import AlleleResult as AlleleResult
from muc_one_span.ladder_clusters import _find_clusters as _find_clusters
from muc_one_span.ladder_clusters import parse_idxstats as parse_idxstats
from muc_one_span.length_candidates import (
    _length_selection_evidence,
    split_cluster_by_read_length,
)
from muc_one_span.read_dominance import (
    evaluate_candidate_pair_dominance,
    extract_read_scores_for_contigs,
)
from muc_one_span.run_status import InsufficientEvidenceError
from muc_one_span.settings import DEFAULT_SETTINGS, AlleleSelectionSettings, ReferenceLayoutSettings
from muc_one_span.tools import run_tool_iter

logger = logging.getLogger(__name__)

# Number of fixed repeat units in the ladder reference (pre-repeats 1-5 + after-repeats 6-9).
# Each contig_N has N canonical X repeats plus these fixed repeats,
# so total allele length = N + PRE_AFTER_REPEAT_COUNT.
PRE_AFTER_REPEAT_COUNT = DEFAULT_SETTINGS.reference_layout.fixed_repeat_count


def _parse_cigar_indel_bp(cigar: str) -> int:
    """Sum of insertion and deletion bases from a CIGAR string."""
    total = 0
    for m in re.finditer(r"(\d+)([ID])", cigar):
        total += int(m.group(1))
    return total


def refine_peak_contig(
    bam_path: Path,
    cluster_contigs: list[str],
    *,
    metric: str = "auto",
    platform: str = "hifi",
) -> dict:
    """Select the best contig from a cluster using alignment quality metrics.

    Scans all reads mapped to the cluster contigs and computes per-contig
    mean alignment score (AS tag) and mean indel length (from CIGAR).
    For HiFi, the contig with the highest mean AS is selected.
    For ONT (or metric='indel'), the contig with the minimum mean indel bp
    among supported contigs is selected to correct for homopolymer drift.
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
                val_str = tag[5:]
                if val_str.lstrip("-").isdigit():
                    contig_stats[contig]["as_sum"] += int(val_str)
                break

    # Compute means and pick best
    metrics: dict[str, dict] = {}
    valid_contigs = [c for c, stats in contig_stats.items() if stats["count"] > 0]
    best_contig = cluster_contigs[0]

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

    use_indel = metric == "indel" or (metric == "auto" and platform == "ont")
    if valid_contigs:
        max_reads = max(contig_stats[c]["count"] for c in valid_contigs)
        threshold = max(3, int(0.25 * max_reads))
        supported = [
            c for c in valid_contigs if contig_stats[c]["count"] >= threshold
        ] or valid_contigs
        if use_indel:
            best_contig = min(
                supported,
                key=lambda c: (
                    contig_stats[c]["indel_sum"] / contig_stats[c]["count"],
                    -contig_stats[c]["as_sum"] / contig_stats[c]["count"],
                ),
            )
        else:
            best_contig = max(
                supported,
                key=lambda c: (
                    contig_stats[c]["as_sum"] / contig_stats[c]["count"],
                    -contig_stats[c]["indel_sum"] / contig_stats[c]["count"],
                ),
            )

    logger.debug("Refined peak contig: %s (metric=%s)", best_contig, "indel" if use_indel else "as")
    return {"best_contig": best_contig, "metrics": metrics}


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

    # Filter valleys: prefer biological candidates (c >= 10) if at least two exist
    candidate_valleys = [v for v in valleys if v[0] >= 10]
    if len(candidate_valleys) < 2:
        candidate_valleys = valleys

    # Take the two lowest valleys normalized by contig reference length
    candidate_valleys.sort(key=lambda x: x[1] / ((x[0] + 9) * 60))
    best_two = sorted(candidate_valleys[:2], key=lambda x: x[0])
    v1, v2 = best_two[0][0], best_two[1][0]

    # Verify the valleys are meaningfully separated (at least 3 contigs apart)
    if abs(v2 - v1) < settings.valley_min_separation:
        return None

    # Split cluster contigs into two sub-clusters by nearest valley
    contigs_dict = dict(cluster["contigs"])
    sub1 = [(c, contigs_dict[c]) for c in sorted(contigs_dict) if abs(c - v1) <= abs(c - v2)]
    sub2 = [(c, contigs_dict[c]) for c in sorted(contigs_dict) if abs(c - v2) < abs(c - v1)]

    if not any(c == v1 for c, _ in sub1):
        sub1.append((v1, 1))
        sub1.sort(key=lambda x: x[0])
    if not any(c == v2 for c, _ in sub2):
        sub2.append((v2, 1))
        sub2.sort(key=lambda x: x[0])

    def _make_sub_cluster(contigs: list[tuple[int, int]], peak_c: int) -> dict:
        total = sum(r for _, r in contigs)
        return {
            "center": peak_c,
            "total_reads": total,
            "contigs": contigs,
            "split_diagnostics": {
                "splitter": "indel_valley",
                "valleys": [v1, v2],
                "target_c": peak_c,
            },
        }

    return [_make_sub_cluster(sub1, v1), _make_sub_cluster(sub2, v2)]


def _build_allele_info(
    cluster: dict,
    best_contig: str | None = None,
    *,
    settings: AlleleSelectionSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
    platform: str = "hifi",
) -> dict:
    """Combine cluster properties, contig names, and repeat conversions.

    Args:
        cluster: Cluster dict with center, total_reads, contigs.
        best_contig: Contig name selected by refine_peak_contig.
            If None, falls back to the weighted center.
        settings: Optional allele selection settings.
        reference_layout: Optional reference layout settings.
        platform: 'hifi' or 'ont'.
    """
    settings = settings or DEFAULT_SETTINGS.allele_selection
    layout = reference_layout or DEFAULT_SETTINGS.reference_layout
    canonical = cluster["center"]
    contig_name = best_contig if best_contig is not None else f"contig_{canonical}"

    # When refine_peak_contig has identified a specific best contig,
    # derive canonical_repeats from its name rather than the cluster
    # center. ONT reads produce wider distributions that can shift
    # the weighted center by ±1-2 vs the indel-refined best contig.
    if best_contig is not None:
        match = re.search(r"_(\d+)$", best_contig)
        if match:
            refined_canonical = int(match.group(1))
            max_shift = (
                max(settings.refinement_max_shift, 2)
                if platform == "ont"
                else settings.refinement_max_shift
            )
            if abs(refined_canonical - canonical) <= max_shift:
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


def detect_alleles(
    counts: dict[int, int],
    min_coverage: int = DEFAULT_SETTINGS.run.min_coverage,
    bam_path: Path | None = None,
    *,
    settings: AlleleSelectionSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
    platform: str = "hifi",
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
        Dictionary with keys ``allele_1``, ``allele_2``, and ``homozygous``, containing
        allele lengths, alignment counts, contig names, and selection evidence.

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
        refined = refine_peak_contig(
            bam_path,
            contig_names,
            metric=settings.refinement_metric,
            platform=platform,
        )
        fit_metrics.update(refined["metrics"])
        best: str = refined["best_contig"]
        return best

    delta_val = settings.score_margin_hifi if platform == "hifi" else settings.score_margin_ont
    min_reads_val = (
        settings.min_dominant_reads_hifi if platform == "hifi" else settings.min_dominant_reads_ont
    )
    min_ratio_val = settings.min_dominance_ratio

    # If only one cluster found but BAM is available, try indel-valley splitting
    if len(clusters) == 1 and bam_path is not None:
        primary_peak_contig = _get_best_contig(clusters[0]) or f"contig_{clusters[0]['center']}"
        split_attempts = [
            split_cluster_by_read_length(
                bam_path, clusters[0], platform=platform, run_tool_iter_func=run_tool_iter
            ),
            _split_cluster_by_indel(bam_path, clusters[0], settings=settings),
        ]
        for sub_clusters in split_attempts:
            if sub_clusters is None:
                continue
            sub_clusters.sort(key=lambda sc: sc["total_reads"], reverse=True)
            sc1, sc2 = sub_clusters[0], sub_clusters[1]
            c1_name = _get_best_contig(sc1) or f"contig_{sc1['center']}"
            c2_name = _get_best_contig(sc2) or f"contig_{sc2['center']}"
            c2_primary = sum(
                fit_metrics.get(f"contig_{c}", {}).get("primary_alignment_records", 0)
                for c, _ in sc2["contigs"]
            )
            if c1_name != c2_name and c2_primary >= min_reads_val:
                sub_scores = extract_read_scores_for_contigs(
                    bam_path, [c1_name, c2_name], run_tool_iter_func=run_tool_iter
                )
                if sub_scores:
                    dom = evaluate_candidate_pair_dominance(
                        sub_scores,
                        c1_name,
                        c2_name,
                        platform=platform,
                        delta=delta_val,
                        min_dominant_reads=min_reads_val,
                        min_ratio=min_ratio_val,
                        c2_primary_records=c2_primary,
                    )
                    if dom.is_valid_second_allele:
                        clusters = sub_clusters
                        break
                    logger.info(
                        "Split candidate rejected by read dominance: %s vs %s (%s)",
                        c1_name,
                        c2_name,
                        dom.rejection_reason,
                    )

        if len(clusters) == 1:
            c1_center = clusters[0]["center"]
            min_gap = settings.min_gap if settings else DEFAULT_SETTINGS.allele_selection.min_gap
            minority_counts = {
                c: r for c, r in int_counts.items() if r >= 3 and abs(c - c1_center) >= min_gap
            }
            if minority_counts:
                minority_sub_clusters = _find_clusters(
                    minority_counts, min_coverage=3, min_gap=min_gap
                )
                for sc in minority_sub_clusters:
                    c1_name = primary_peak_contig
                    c2_name = _get_best_contig(sc) or f"contig_{sc['center']}"
                    c2_primary = sum(
                        fit_metrics.get(f"contig_{c}", {}).get("primary_alignment_records", 0)
                        for c, _ in sc["contigs"]
                    )
                    if c2_primary < min_reads_val:
                        logger.info(
                            "Minority candidate %s rejected: insufficient primary records (%d < %d)",
                            c2_name,
                            c2_primary,
                            min_reads_val,
                        )
                        continue
                    sub_scores = extract_read_scores_for_contigs(
                        bam_path,
                        [c1_name, c2_name],
                        run_tool_iter_func=run_tool_iter,
                    )
                    if sub_scores:
                        dom_sub = evaluate_candidate_pair_dominance(
                            sub_scores,
                            c1_name,
                            c2_name,
                            platform=platform,
                            delta=delta_val,
                            min_dominant_reads=min_reads_val,
                            min_ratio=min_ratio_val,
                            c2_primary_records=c2_primary,
                        )
                        if dom_sub.is_valid_second_allele:
                            clusters.append(sc)
                            break

    selection_evidence = _length_selection_evidence(
        int_counts,
        min_coverage,
        bam_path,
        clusters[2:],
        run_tool_iter_func=run_tool_iter,
    )

    def _with_support(cluster: dict) -> dict:
        info = _build_allele_info(
            cluster,
            _get_best_contig(cluster),
            settings=settings,
            reference_layout=reference_layout,
            platform=platform,
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

    allele_2 = _with_support(clusters[1]) if len(clusters) >= 2 else None

    if allele_2 is None or allele_1["length"] == allele_2["length"]:
        if allele_2 is not None:
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
