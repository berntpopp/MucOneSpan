# tests/unit/test_read_dominance.py
"""Tests for read dominance scoring and candidate validation."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from muc_one_span.read_dominance import (
    evaluate_candidate_pair_dominance,
    extract_read_scores_for_contigs,
)


class TestEvaluateCandidatePairDominance:
    """Tests for pairwise read dominance evaluation."""

    def test_homozygous_rejection_zero_d2(self):
        """When all reads prefer c1, c2 is rejected."""
        read_scores = {
            f"read_{i}": {"contig_51": 1600, "contig_54": 1500}  # diff = -100 (< -43)
            for i in range(50)
        }
        res = evaluate_candidate_pair_dominance(read_scores, "contig_51", "contig_54", "hifi")
        assert res.is_valid_second_allele is False
        assert res.d1_count == 50
        assert res.d2_count == 0
        assert res.ambiguous_count == 0
        assert "Insufficient dominant reads" in str(res.rejection_reason)

    def test_asymmetric_minority_allele_accepted(self):
        """4 dominant minority reads out of 54 are accepted on HiFi (floor=3)."""
        read_scores = {f"maj_read_{i}": {"contig_16": 1200} for i in range(50)}
        for i in range(4):
            # minority reads prefer contig_131 by huge margin
            read_scores[f"min_read_{i}"] = {"contig_16": 200, "contig_131": 1800}

        res = evaluate_candidate_pair_dominance(
            read_scores, "contig_16", "contig_131", "hifi", min_ratio=0.03
        )
        assert res.is_valid_second_allele is True
        assert res.d1_count == 50
        assert res.d2_count == 4
        assert res.is_valid_second_allele is True

    def test_low_support_below_floor_rejected(self):
        """Only 2 reads for c2 is below min_dominant_reads (3 for HiFi)."""
        read_scores = {f"maj_read_{i}": {"contig_51": 1600, "contig_52": 1500} for i in range(40)}
        read_scores["min_read_1"] = {"contig_51": 1500, "contig_52": 1600}
        read_scores["min_read_2"] = {"contig_51": 1500, "contig_52": 1600}

        res = evaluate_candidate_pair_dominance(read_scores, "contig_51", "contig_52", "hifi")
        assert res.is_valid_second_allele is False
        assert res.d2_count == 2

    def test_ambiguous_reads_not_counted_as_dominant(self):
        """Reads with score difference <= delta are counted as ambiguous."""
        read_scores = {
            "read_1": {"c1": 1000, "c2": 1020},  # diff = 20 <= 43 -> ambiguous
            "read_2": {"c1": 1000, "c2": 990},  # diff = -10 >= -43 -> ambiguous
            "read_3": {"c1": 1000, "c2": 1100},  # diff = 100 > 43 -> D2
            "read_4": {"c1": 1000, "c2": 1100},  # diff = 100 > 43 -> D2
            "read_5": {"c1": 1000, "c2": 1100},  # diff = 100 > 43 -> D2
        }
        res = evaluate_candidate_pair_dominance(read_scores, "c1", "c2", "hifi")
        assert res.ambiguous_count == 2
        assert res.d2_count == 3
        assert res.d1_count == 0
        assert res.is_valid_second_allele is True

    def test_close_candidate_missing_score_is_ambiguous(self):
        """Close candidate contigs (<10 repeats apart) without score on c1 remain ambiguous."""
        # contig_71 (80 repeats) vs contig_66 (75 repeats): delta = 5 repeats
        read_scores = {f"maj_{i}": {"contig_71": 4200, "contig_66": 3900} for i in range(50)}
        # 3 reads only aligned to contig_66 (missing contig_71 score)
        read_scores["tail_1"] = {"contig_66": 4150}
        read_scores["tail_2"] = {"contig_66": 4100}
        read_scores["tail_3"] = {"contig_66": 3950}

        res = evaluate_candidate_pair_dominance(read_scores, "contig_71", "contig_66", "hifi")
        # Tail reads must be ambiguous, not d2!
        assert res.d2_count == 0
        assert res.ambiguous_count == 3
        assert res.is_valid_second_allele is False

    def test_extreme_asymmetry_mathematical_bound_accepts_dominant(self):
        """Distant candidate contig (>=10 repeats apart) with score exceeding theoretical max is decisive."""
        # contig_16 (25 repeats, length 2500, max score 5000) vs contig_131 (140 repeats)
        read_scores = {f"maj_{i}": {"contig_16": 1300} for i in range(50)}
        for i in range(4):
            # score 7500 > 5000 + 43
            read_scores[f"min_{i}"] = {"contig_131": 7500}

        res = evaluate_candidate_pair_dominance(read_scores, "contig_16", "contig_131", "hifi")
        assert res.d2_count == 4
        assert res.is_valid_second_allele is True


class TestExtractReadScoresForContigs:
    """Tests for extracting alignment scores from BAM via samtools view."""

    def test_extract_scores_parses_as_tag(self, tmp_path: Path):
        mock_sam = [
            "read1\t0\tcontig_51\t100\t60\t60M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:1650",
            "read1\t256\tcontig_54\t100\t0\t60M\t*\t0\t0\tACGT\tIIII\tNM:i:3\tAS:i:1520",
            "read2\t0\tcontig_51\t100\t60\t60M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:1640",
        ]
        bam = tmp_path / "mock.bam"
        bam.touch()

        with patch("muc_one_span.read_dominance.run_tool_iter", return_value=iter(mock_sam)):
            scores = extract_read_scores_for_contigs(bam, ["contig_51", "contig_54"])

        assert "read1" in scores
        assert scores["read1"]["contig_51"] == 1650
        assert scores["read1"]["contig_54"] == 1520
        assert scores["read2"]["contig_51"] == 1640
        assert "contig_54" not in scores["read2"]
