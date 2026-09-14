# tests/unit/test_ladder.py
"""Tests for reference ladder generation."""

from __future__ import annotations

import pytest

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.ladder import build_contig, generate_ladder_fasta


@pytest.fixture
def repeat_dict():
    return load_repeat_dictionary()


class TestBuildContig:
    """Tests for building a single ladder contig."""

    def test_contig_contains_flanking(self, repeat_dict):
        """Contig starts with left flank and ends with right flank."""
        contig = build_contig(5, repeat_dict, flank_length=100)
        assert contig["sequence"].startswith(repeat_dict.flanking_left[:100])
        assert contig["sequence"].endswith(repeat_dict.flanking_right[:100])

    def test_contig_contains_pre_and_after_repeats(self, repeat_dict):
        """Contig contains pre-repeats and after-repeats."""
        contig = build_contig(1, repeat_dict, flank_length=0)
        seq = contig["sequence"]
        # Pre-repeats 1,2,3,4,5 should be at the start
        pre_concat = "".join(repeat_dict.repeats[rid] for rid in ["1", "2", "3", "4", "5"])
        assert seq.startswith(pre_concat)
        # After-repeats 6,7,8,9 should be at the end
        after_concat = "".join(repeat_dict.repeats[rid] for rid in ["6", "7", "8", "9"])
        assert seq.endswith(after_concat)

    def test_contig_repeat_count(self, repeat_dict):
        """Contig with N=10 has 10 canonical X repeats in the middle."""
        contig = build_contig(10, repeat_dict, flank_length=0)
        x_seq = repeat_dict.repeats["X"]
        pre_len = sum(len(repeat_dict.repeats[rid]) for rid in ["1", "2", "3", "4", "5"])
        after_len = sum(len(repeat_dict.repeats[rid]) for rid in ["6", "7", "8", "9"])
        expected_len = pre_len + (10 * len(x_seq)) + after_len
        assert len(contig["sequence"]) == expected_len

    def test_contig_name(self, repeat_dict):
        """Contig name follows expected pattern."""
        contig = build_contig(42, repeat_dict, flank_length=100)
        assert contig["name"] == "contig_42"


class TestGenerateLadderFasta:
    """Tests for generating the full ladder FASTA."""

    def test_generates_correct_count(self, repeat_dict, tmp_path):
        """Ladder with range 1-5 produces 5 contigs."""
        output = tmp_path / "ladder.fa"
        generate_ladder_fasta(repeat_dict, output, min_units=1, max_units=5, flank_length=100)
        content = output.read_text()
        headers = [line for line in content.splitlines() if line.startswith(">")]
        assert len(headers) == 5

    def test_default_range_generates_150_contigs(self, repeat_dict, tmp_path):
        """Default range 1-150 produces 150 contigs."""
        output = tmp_path / "ladder.fa"
        generate_ladder_fasta(repeat_dict, output, flank_length=100)
        content = output.read_text()
        headers = [line for line in content.splitlines() if line.startswith(">")]
        assert len(headers) == 150

    def test_fasta_format(self, repeat_dict, tmp_path):
        """Output is valid FASTA format."""
        output = tmp_path / "ladder.fa"
        generate_ladder_fasta(repeat_dict, output, min_units=1, max_units=2, flank_length=50)
        content = output.read_text()
        lines = content.strip().splitlines()
        assert lines[0].startswith(">contig_1")
        # Sequence lines should only contain ACGT
        for line in lines:
            if not line.startswith(">"):
                assert all(c in "ACGTacgt" for c in line), f"Invalid chars in: {line[:50]}"
