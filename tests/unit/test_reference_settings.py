"""Configured reference layout uses selected IDs, preserving explicit overrides."""

from dataclasses import replace
from pathlib import Path

import pytest

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.ladder import build_contig, generate_ladder_fasta
from muc_one_span.settings import ConsensusSettings, ReferenceLayoutSettings


def test_layout_uses_selected_ids_in_order_not_dictionary_categories() -> None:
    rd = load_repeat_dictionary()
    layout = ReferenceLayoutSettings(pre=("5", "1"), after=("9",))
    result = build_contig(2, rd, 0, reference_layout=layout)
    assert result["sequence"] == "".join(rd.repeats[rid] for rid in ("5", "1", "X", "X", "9"))
    assert len(result["sequence"]) == (2 + layout.fixed_repeat_count) * rd.repeat_length_bp


def test_missing_selected_repeat_is_rejected_before_output(tmp_path: Path) -> None:
    rd = load_repeat_dictionary()
    rd = replace(rd, repeats={key: value for key, value in rd.repeats.items() if key != "1"})
    output = tmp_path / "missing.fa"
    with pytest.raises(ValueError, match="1"):
        generate_ladder_fasta(rd, output)
    assert not output.exists()


def test_flank_setting_and_explicit_override(tmp_path: Path) -> None:
    rd = load_repeat_dictionary()
    settings = ConsensusSettings(flank_length=7)
    configured = build_contig(1, rd, settings=settings)["sequence"]
    assert configured == build_contig(1, rd, 7)["sequence"]
    assert build_contig(1, rd, 0, settings=settings) == build_contig(1, rd, 0)
    output = generate_ladder_fasta(rd, tmp_path / "configured.fa", 1, 1, settings=settings)
    assert "".join(output.read_text().splitlines()[1:]) == configured


def test_ladder_range_setting_and_explicit_override(tmp_path: Path) -> None:
    rd = load_repeat_dictionary()
    settings = ReferenceLayoutSettings(min_units=2, max_units=3)
    output = generate_ladder_fasta(rd, tmp_path / "configured.fa", reference_layout=settings)
    assert [line for line in output.read_text().splitlines() if line.startswith(">")] == [
        ">contig_2",
        ">contig_3",
    ]
    output = generate_ladder_fasta(rd, tmp_path / "explicit.fa", 1, 1, reference_layout=settings)
    assert [line for line in output.read_text().splitlines() if line.startswith(">")] == [
        ">contig_1"
    ]


@pytest.mark.parametrize("left,right", [("AAAA", "GGGG"), ("AAAA", "GGGGGG"), ("AAAAAA", "GGGG")])
def test_requested_flanks_must_fit_both_dictionary_sequences(
    tmp_path: Path, left: str, right: str
) -> None:
    rd = replace(load_repeat_dictionary(), flanking_left=left, flanking_right=right)
    settings = ConsensusSettings(flank_length=6)
    with pytest.raises(ValueError, match=r"flank_length.*available"):
        build_contig(1, rd, settings=settings)
    output = tmp_path / "invalid.fa"
    with pytest.raises(ValueError, match=r"flank_length.*available"):
        generate_ladder_fasta(rd, output, settings=settings)
    assert not output.exists()


def test_explicit_valid_flank_override_precedes_oversized_setting() -> None:
    rd = replace(load_repeat_dictionary(), flanking_left="AAAA", flanking_right="GGGG")
    settings = ConsensusSettings(flank_length=6)
    assert build_contig(1, rd, 4, settings=settings) == build_contig(1, rd, 4)
    empty = replace(rd, flanking_left="", flanking_right="")
    assert build_contig(1, empty, 0, settings=settings) == build_contig(1, rd, 0)
