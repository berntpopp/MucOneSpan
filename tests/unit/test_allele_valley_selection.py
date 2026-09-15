"""Tests for normalized indel-rate valley selection and fragment filtering."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from muc_one_span.alleles import _split_cluster_by_indel
from muc_one_span.settings import AlleleSelectionSettings


def _mock_sam_line(name: str, contig: int, cigar: str) -> str:
    return f"{name}\t0\tcontig_{contig}\t1\t60\t{cigar}\t*\t0\t0\tAAAA\tIIII\tAS:i:60"


def test_valley_selection_prefers_biological_alleles_over_short_fragments(tmp_path: Path) -> None:
    """Valleys at c=41 and c=69 must be selected over c=5 fragment even if c=5 has lower raw indel bp."""
    # Build a cluster containing contigs from 1 to 75
    contigs = [(c, 100) for c in range(1, 75)]
    cluster = {"center": 45, "total_reads": len(contigs) * 100, "contigs": contigs}

    # Simulate SAM records:
    # Contig 5: mean indel = 75 bp (ref length ~14*60 = 840 bp, rate ~0.089)
    # Contig 41: mean indel = 72 bp (ref length ~50*60 = 3000 bp, rate ~0.024)
    # Contig 69: mean indel = 110 bp (ref length ~78*60 = 4680 bp, rate ~0.0235)
    # Other contigs have higher indel values creating peaks around these valleys
    sam_lines: list[str] = []
    for c in range(1, 75):
        if c == 5:
            cigar = "60M75I"
        elif c == 41:
            cigar = "60M72I"
        elif c == 69:
            cigar = "60M110I"
        elif c < 10:
            cigar = "60M120I"
        elif c < 50:
            cigar = "60M180I"
        else:
            cigar = "60M200I"
        sam_lines.append(_mock_sam_line(f"read_{c}", c, cigar))

    with patch("muc_one_span.alleles.run_tool_iter", return_value=iter(sam_lines)):
        sub_clusters = _split_cluster_by_indel(
            tmp_path / "mock.bam",
            cluster,
            settings=AlleleSelectionSettings(valley_min_points=3, valley_min_separation=3),
        )

    assert sub_clusters is not None, "Expected cluster to be split into two sub-clusters"
    valleys = sub_clusters[0]["split_diagnostics"]["valleys"]
    assert valleys == [41, 69], f"Expected valleys [41, 69], but got {valleys}"


def test_valley_selection_handles_cases_with_all_small_contigs(tmp_path: Path) -> None:
    """When only c < 10 contigs exist, valleys are not dropped."""
    contigs = [(c, 50) for c in range(1, 9)]
    cluster = {"center": 4, "total_reads": len(contigs) * 50, "contigs": contigs}

    sam_lines: list[str] = []
    for c in range(1, 9):
        if c == 2:
            cigar = "60M20I"
        elif c == 6:
            cigar = "60M25I"
        else:
            cigar = "60M80I"
        sam_lines.append(_mock_sam_line(f"read_{c}", c, cigar))

    with patch("muc_one_span.alleles.run_tool_iter", return_value=iter(sam_lines)):
        sub_clusters = _split_cluster_by_indel(
            tmp_path / "mock.bam",
            cluster,
            settings=AlleleSelectionSettings(valley_min_points=3, valley_min_separation=3),
        )

    assert sub_clusters is not None
    valleys = sub_clusters[0]["split_diagnostics"]["valleys"]
    assert valleys == [2, 6]
