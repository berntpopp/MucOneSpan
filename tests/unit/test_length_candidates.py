# tests/unit/test_length_candidates.py
"""Tests for length candidate discovery and validation."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from muc_one_span.length_candidates import (
    CandidateCluster,
    _length_selection_evidence,
    discover_length_candidates,
    validate_candidates_with_dominance,
)


def test_discover_length_candidates_empty() -> None:
    assert discover_length_candidates({}) == []
    assert discover_length_candidates({50: 2}, discovery_floor=3) == []


def test_discover_length_candidates_single_and_secondary() -> None:
    counts = {51: 40, 52: 2, 71: 15}
    candidates = discover_length_candidates(counts, min_coverage=10)
    assert len(candidates) == 2
    assert candidates[0].canonical_repeats == 51
    assert candidates[0].total_reads == 40
    assert candidates[1].canonical_repeats == 71
    assert candidates[1].total_reads == 15


def test_discover_length_candidates_minority_with_primary() -> None:
    counts = {51: 50, 131: 4}
    prim = {51: 50, 131: 3}
    candidates = discover_length_candidates(
        counts, primary_counts=prim, min_coverage=10, discovery_floor=3
    )
    assert len(candidates) == 2
    assert candidates[1].canonical_repeats == 131
    assert candidates[1].primary_reads == 3


def test_validate_candidates_with_dominance_edge_cases(tmp_path: Path) -> None:
    c1 = CandidateCluster(51, "contig_51", 40, 40, ["contig_51"])
    # Same candidate
    valid, _ = validate_candidates_with_dominance(tmp_path / "x.bam", c1, c1)
    assert not valid

    # None candidate_2
    valid, _ = validate_candidates_with_dominance(tmp_path / "x.bam", c1, None)
    assert not valid

    # BAM does not exist
    c2 = CandidateCluster(71, "contig_71", 15, 15, ["contig_71"])
    valid, _ = validate_candidates_with_dominance(tmp_path / "nonexistent.bam", c1, c2)
    assert valid

    # BAM with no scores
    bam = tmp_path / "empty.bam"
    bam.touch()
    with patch("muc_one_span.length_candidates.extract_read_scores_for_contigs", return_value={}):
        valid, _ = validate_candidates_with_dominance(bam, c1, c2)
        assert not valid


def test_length_selection_evidence_without_bam() -> None:
    counts = {51: 30, 71: 4}
    res = _length_selection_evidence(counts, min_coverage=10, bam_path=None, unselected_clusters=[])
    assert res["minimum_coverage"] == 10
    assert res["excluded_subthreshold_contigs"][0]["contig_name"] == "contig_71"
    assert res["excluded_subthreshold_primary_alignment_records"] is None


def test_split_cluster_by_read_length_delta_2(tmp_path: Path) -> None:
    """Verify split_cluster_by_read_length correctly splits Delta=2 (e.g. 40/42) without midpoint shift."""
    from muc_one_span.length_candidates import split_cluster_by_read_length

    # 40 repeats: 40 * 60 + 30 = 2430 bp (c1 = 40 - 9 = 31)
    # 42 repeats: 42 * 60 + 30 = 2550 bp (c2 = 42 - 9 = 33)
    # Simulated SAM lines from runner
    # 20 reads of length 2430 bp, 20 reads of length 2550 bp
    seq_40 = "A" * 2430
    seq_42 = "A" * 2550
    lines = []
    for i in range(20):
        lines.append(f"read_40_{i}\t0\tcontig_31\t1\t60\t2430M\t*\t0\t0\t{seq_40}\t*")
    for i in range(20):
        lines.append(f"read_42_{i}\t0\tcontig_33\t1\t60\t2550M\t*\t0\t0\t{seq_42}\t*")

    bam = tmp_path / "mock.bam"
    bam.touch()

    cluster = {
        "center": 32,
        "total_reads": 40,
        "contigs": [(30, 5), (31, 15), (32, 10), (33, 10)],
    }

    sub_clusters = split_cluster_by_read_length(
        bam,
        cluster,
        platform="hifi",
        min_reads=5,
        run_tool_iter_func=lambda cmd: lines,
    )

    assert sub_clusters is not None
    assert len(sub_clusters) == 2
    # Centers must be locked to 31 and 33, intermediate contig 32 excluded from determining centers
    assert sub_clusters[0]["center"] == 31
    assert sub_clusters[1]["center"] == 33
    assert sub_clusters[0]["split_diagnostics"]["splitter"] == "read_length"
    assert sub_clusters[0]["split_diagnostics"]["delta"] == 120


def test_split_cluster_by_read_length_monomodal(tmp_path: Path) -> None:
    """Monomodal distribution returns None (no split)."""
    from muc_one_span.length_candidates import split_cluster_by_read_length

    seq_60 = "A" * 3630
    lines = [f"read_{i}\t0\tcontig_51\t1\t60\t3630M\t*\t0\t0\t{seq_60}\t*" for i in range(30)]

    bam = tmp_path / "mock.bam"
    bam.touch()

    cluster = {"center": 51, "total_reads": 30, "contigs": [(51, 30)]}
    res = split_cluster_by_read_length(
        bam, cluster, platform="hifi", min_reads=5, run_tool_iter_func=lambda cmd: lines
    )
    assert res is None
