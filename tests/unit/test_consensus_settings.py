"""Consensus settings control flanks and exact anchors without changing defaults."""

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.consensus import build_consensus_per_allele, trim_flanking
from muc_one_span.settings import ConsensusSettings, ReferenceLayoutSettings


def test_flank_setting_and_explicit_override(tmp_path: Path) -> None:
    source = tmp_path / "full.fa"
    source.write_text(">c\nAACCCGG\n")
    output = tmp_path / "trim.fa"
    settings = ConsensusSettings(flank_length=2)
    trim_flanking(source, None, output, settings=settings)
    assert output.read_text() == ">c_vntr\nCCC\n"
    trim_flanking(source, 1, output, settings=settings)
    assert output.read_text() == ">c_vntr\nACCCG\n"


def test_configured_anchor_ids_components_and_tolerance(tmp_path: Path) -> None:
    rd = replace(
        load_repeat_dictionary(),
        repeats={"A": "CGTT", "B": "GAAC"},
        flanking_left="TTAA",
        flanking_right="CCGG",
    )
    layout = ReferenceLayoutSettings(pre=("A",), after=("B",))
    source = tmp_path / "full.fa"
    source.write_text(">c\nTTTAACGTTGAACCCGG\n")
    output = tmp_path / "trim.fa"
    context = {}
    trim_flanking(
        source,
        None,
        output,
        rd,
        context=context,
        settings=ConsensusSettings(flank_length=4, anchor_bases=2, anchor_tolerance=2),
        reference_layout=layout,
    )
    assert output.read_text() == ">c_vntr\nCGTTGAAC\n"
    assert context["left_trim_method"] == context["right_trim_method"] == "exact_anchor"
    trim_flanking(
        source,
        None,
        output,
        rd,
        settings=ConsensusSettings(flank_length=4, anchor_bases=2, anchor_tolerance=0),
        reference_layout=layout,
    )
    assert output.read_text() == ">c_vntr\nACGTTGAAC\n"


def test_selected_anchor_ids_are_validated_even_without_flanks(tmp_path: Path) -> None:
    source = tmp_path / "full.fa"
    source.write_text(">c\nACGT\n")
    with pytest.raises(ValueError, match="missing"):
        trim_flanking(
            source,
            0,
            tmp_path / "out.fa",
            load_repeat_dictionary(),
            reference_layout=ReferenceLayoutSettings(pre=("missing",), after=("9",)),
        )


def test_per_allele_pipeline_passes_flank_settings(tmp_path: Path) -> None:
    def consensus(_reference, _vcf, output, **_kwargs):
        output.write_text(">contig_1\nAACCCGG\n")
        return output

    rd = replace(
        load_repeat_dictionary(),
        repeats={"A": "C", "B": "C"},
        flanking_left="AA",
        flanking_right="GG",
    )
    layout = ReferenceLayoutSettings(pre=("A",), after=("B",))
    alleles = {"allele_1": {"contig_name": "contig_1", "length": 10}}
    with (
        patch("muc_one_span.consensus.run_tool", return_value=">contig_1\nAACCCGG\n"),
        patch("muc_one_span.consensus.select_vcf_sample", return_value="sample"),
        patch("muc_one_span.consensus.build_consensus", side_effect=consensus),
    ):
        outputs = build_consensus_per_allele(
            tmp_path / "ref.fa",
            {"allele_1": tmp_path / "v.vcf"},
            alleles,
            tmp_path / "out",
            settings=ConsensusSettings(flank_length=2, anchor_bases=2, anchor_tolerance=0),
            repeat_dict=rd,
            reference_layout=layout,
        )
    assert outputs["allele_1"].read_text() == ">contig_1_vntr\nCCC\n"


@pytest.mark.parametrize("anchor_bases", [2, 20, 100])
def test_sliced_anchor_components_preserve_short_repeat_boundaries(
    tmp_path: Path, anchor_bases: int
) -> None:
    rd = replace(
        load_repeat_dictionary(),
        repeats={"A": "CGTT", "B": "GAAC"},
        flanking_left="TTAA",
        flanking_right="CCGG",
    )
    source = tmp_path / "full.fa"
    source.write_text(">c\nTTAACGTTGAACCCGG\n")
    output = tmp_path / "trim.fa"
    context = {}
    trim_flanking(
        source,
        4,
        output,
        rd,
        context=context,
        settings=ConsensusSettings(anchor_bases=anchor_bases, anchor_tolerance=100),
        reference_layout=ReferenceLayoutSettings(pre=("A",), after=("B",)),
    )
    assert output.read_text() == ">c_vntr\nCGTTGAAC\n"
    assert context["trim_start"] == 4
    assert context["trim_end"] == 12
    assert context["left_trim_method"] == context["right_trim_method"] == "exact_anchor"


def test_zero_tolerance_finds_anchors_at_expected_boundaries(tmp_path: Path) -> None:
    rd = replace(
        load_repeat_dictionary(),
        repeats={"A": "CGTT", "B": "GAAC"},
        flanking_left="TTAA",
        flanking_right="CCGG",
    )
    source = tmp_path / "full.fa"
    source.write_text(">c\nTTAACGTTGAACCCGG\n")
    context = {}
    trim_flanking(
        source,
        4,
        tmp_path / "trim.fa",
        rd,
        context=context,
        settings=ConsensusSettings(anchor_bases=2, anchor_tolerance=0),
        reference_layout=ReferenceLayoutSettings(pre=("A",), after=("B",)),
    )
    assert context == {
        "trim_start": 4,
        "trim_end": 12,
        "left_trim_method": "exact_anchor",
        "right_trim_method": "exact_anchor",
    }


def test_short_flank_anchor_uses_same_prefix_as_ladder(tmp_path: Path) -> None:
    rd = replace(
        load_repeat_dictionary(),
        repeats={"A": "CGTT", "B": "GAAC"},
        flanking_left="TTAA",
        flanking_right="CCGG",
    )
    source = tmp_path / "full.fa"
    # A one-base insertion in the left prefix retained by the ladder.
    source.write_text(">c\nTTTCGTTGAACCC\n")
    output = tmp_path / "trim.fa"
    context = {}
    trim_flanking(
        source,
        2,
        output,
        rd,
        context=context,
        settings=ConsensusSettings(anchor_bases=2, anchor_tolerance=1, proximal_flank=False),
        reference_layout=ReferenceLayoutSettings(pre=("A",), after=("B",)),
    )
    assert output.read_text() == ">c_vntr\nCGTTGAAC\n"
    assert context["trim_start"] == 3
    assert context["trim_end"] == 11


def test_default_prefix_anchor_corrects_one_base_flank_insertion(tmp_path: Path) -> None:
    rd = load_repeat_dictionary()
    vntr = rd.repeats["1"] + rd.repeats["X"] + rd.repeats["9"]
    left = rd.flanking_left[-500:]
    source = tmp_path / "full.fa"
    source.write_text(
        ">c\n" + left[:250] + "A" + left[250:] + vntr + rd.flanking_right[:500] + "\n"
    )
    output = tmp_path / "trim.fa"
    context = {}
    trim_flanking(source, None, output, rd, context=context)
    assert output.read_text() == ">c_vntr\n" + vntr + "\n"
    assert context["trim_start"] == 501
    assert context["trim_end"] == 681


@pytest.mark.parametrize("side", ["left", "right"])
@pytest.mark.parametrize("shift", [-1, 1])
@pytest.mark.parametrize("tolerance", [0, 1])
def test_tolerance_bounds_boundary_displacement_on_both_sides(
    tmp_path: Path, side: str, shift: int, tolerance: int
) -> None:
    rd = replace(
        load_repeat_dictionary(),
        repeats={"A": "GTTG", "B": "CACA"},
        flanking_left="AAAA",
        flanking_right="GGGG",
    )
    left, right = rd.flanking_left, rd.flanking_right
    if side == "left":
        left = "T" + left if shift == 1 else left[1:]
    else:
        right = right + "T" if shift == 1 else right[:-1]
    source = tmp_path / "full.fa"
    source.write_text(">c\n" + left + "GTTGCACA" + right + "\n")
    output = tmp_path / "trim.fa"
    context = {}
    trim_flanking(
        source,
        4,
        output,
        rd,
        context=context,
        settings=ConsensusSettings(anchor_bases=2, anchor_tolerance=tolerance),
        reference_layout=ReferenceLayoutSettings(pre=("A",), after=("B",)),
    )
    assert context[side + "_trim_method"] == (
        "exact_anchor" if tolerance else "fixed_anchor_not_found"
    )
    if tolerance:
        assert output.read_text() == ">c_vntr\nGTTGCACA\n"


def test_oversized_flank_extent_rejects_silent_vntr_truncation(tmp_path: Path) -> None:
    rd = replace(
        load_repeat_dictionary(),
        repeats={"A": "GTTG", "X": "ACGT", "B": "CACA"},
        flanking_left="AAAA",
        flanking_right="GGGG",
    )
    layout = ReferenceLayoutSettings(pre=("A",), after=("B",))
    settings = ConsensusSettings(flank_length=6, anchor_bases=2, anchor_tolerance=0)
    full = tmp_path / "full.fa"
    full.write_text(">c\nAAAAGTTGACGTCACAGGGG\n")
    output = tmp_path / "trim.fa"
    with pytest.raises(ValueError, match=r"flank_length.*available"):
        trim_flanking(full, None, output, rd, settings=settings, reference_layout=layout)
    assert not output.exists()
    trim_flanking(full, 4, output, rd, settings=settings, reference_layout=layout)
    assert output.read_text() == ">c_vntr\nGTTGACGTCACA\n"


def test_invalid_flank_extent_fails_before_consensus_tools(tmp_path: Path) -> None:
    rd = replace(load_repeat_dictionary(), flanking_left="AAAA", flanking_right="GGGG")
    with (
        patch("muc_one_span.consensus.run_tool", side_effect=AssertionError("tool reached")),
        pytest.raises(ValueError, match=r"flank_length.*available"),
    ):
        build_consensus_per_allele(
            tmp_path / "ref.fa",
            {"allele_1": tmp_path / "v.vcf"},
            {"allele_1": {"length": 10}},
            tmp_path / "out",
            repeat_dict=rd,
            settings=ConsensusSettings(flank_length=6),
        )
    assert not (tmp_path / "out").exists()


def test_rebuilt_unresolved_candidate_clears_stale_vntr_phase_status(tmp_path: Path) -> None:
    alleles = {
        "allele_1": {
            "length": 10,
            "contig_name": "contig_1",
            "independent_haplotype_evidence": False,
            "vntr_phase_status": "distinct_genotype_candidates",
        }
    }
    with (
        patch("muc_one_span.consensus.run_tool", return_value=">c\nACGT\n"),
        patch("muc_one_span.consensus.select_vcf_sample", return_value="sample"),
    ):
        result = build_consensus_per_allele(
            tmp_path / "ref.fa",
            {"allele_1": tmp_path / "v.vcf"},
            alleles,
            tmp_path / "out",
            flank_length=0,
        )
    assert result["allele_1"].read_text() == ">c_vntr\nACGT\n"
    assert alleles["allele_1"]["vntr_phase_status"] == "unresolved"
    assert alleles["allele_1"]["independent_haplotype_evidence"] is False
