"""Small actual-reference fixtures shared by classification unit tests."""

from pathlib import Path

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.config import load_repeat_dictionary


@pytest.fixture
def indel_concordance(tmp_path: Path) -> tuple:
    """Three X units with a known two-base insertion in the middle unit."""
    rd = load_repeat_dictionary()
    x = rd.repeats["X"]
    seq = x + x[:40] + "TT" + x[40:] + x
    reference = "A" * 500 + x * 3 + "T" * 500
    full = "A" * 500 + seq + "T" * 500
    ref_path, full_path = tmp_path / "ref.fa", tmp_path / "full.fa"
    ref_path.write_text(">contig_1\n" + reference + "\n")
    full_path.write_text(">contig_1\n" + full + "\n")
    pos = 600
    variants = [
        {
            "chrom": "contig_1",
            "pos": pos,
            "ref": reference[pos - 1],
            "alt": reference[pos - 1] + "TT",
            "qual": 25.0,
            "genotype": "1/1",
        }
    ]
    context = {
        "reference_path": str(ref_path),
        "full_consensus_path": str(full_path),
        "trim_start": 500,
        "trim_end": len(full) - 500,
        "chrom": "contig_1",
        "haplotype": "I",
    }
    return classify_sequence(seq, rd), variants, seq, rd, context
