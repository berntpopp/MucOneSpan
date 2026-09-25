"""Each ladder allele-selection setting (#74) is honoured at a non-default value."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.alleles import (
    _build_allele_info,
    _split_cluster_by_indel,
    detect_alleles,
    refine_peak_contig,
)
from muc_one_span.length_candidates import split_cluster_by_read_length
from muc_one_span.read_dominance import (
    DEFAULT_DELTA_AS,
    DEFAULT_MIN_DOMINANT_READS,
    DEFAULT_MIN_RATIO,
    DominanceScore,
    evaluate_candidate_pair_dominance,
)
from muc_one_span.settings import (
    DEFAULT_SETTINGS,
    AlleleSelectionSettings,
    ReferenceLayoutSettings,
)

DEFAULT = AlleleSelectionSettings()


def _sam(name: str, contig: str, cigar: str = "60M", score: int = 100, seq: str = "A") -> str:
    return f"{name}\t0\t{contig}\t1\t60\t{cigar}\t*\t0\t0\t{seq}\t*\tAS:i:{score}"


# refine_peak_contig: supported-contig threshold -------------------------------------------


@pytest.mark.parametrize(
    ("settings", "best"),
    [
        (None, "contig_40"),
        (
            AlleleSelectionSettings(
                refinement_min_supported_records=1, refinement_supported_fraction=0.01
            ),
            "contig_41",
        ),
    ],
)
def test_refine_peak_contig_support_threshold_is_configurable(
    tmp_path: Path, settings: AlleleSelectionSettings | None, best: str
) -> None:
    lines = [_sam(f"a{i}", "contig_40", score=100) for i in range(40)]
    lines += [_sam(f"b{i}", "contig_41", score=500) for i in range(2)]
    with patch("muc_one_span.alleles.run_tool_iter", return_value=iter(lines)):
        result = refine_peak_contig(
            tmp_path / "m.bam", ["contig_40", "contig_41"], settings=settings
        )
    assert result["best_contig"] == best


# _split_cluster_by_indel: biological-valley floor and layout-normalised ranking ------------


def _indel_cluster(valleys: dict[int, int], last: int) -> tuple[dict, list[str]]:
    contigs = [(c, 10) for c in range(1, last + 1)]
    lines = [_sam(f"r{c}", f"contig_{c}", cigar=f"60M{valleys.get(c, 200)}I") for c, _ in contigs]
    return {"center": last // 2, "total_reads": 10 * last, "contigs": contigs}, lines


def _valley_pair(tmp_path: Path, valleys: dict[int, int], last: int, **kwargs: Any) -> list:
    cluster, lines = _indel_cluster(valleys, last)
    with patch("muc_one_span.alleles.run_tool_iter", return_value=iter(lines)):
        result = _split_cluster_by_indel(tmp_path / "m.bam", cluster, **kwargs)
    assert result is not None
    return list(result[0]["split_diagnostics"]["valleys"])


def test_valley_floor_is_configurable_with_three_valleys(tmp_path: Path) -> None:
    """Contig 6 has the lowest length-normalised indel rate but is below the default floor."""
    valleys = {6: 10, 12: 40, 20: 60}
    assert _valley_pair(tmp_path, valleys, 24) == [12, 20]
    settings = AlleleSelectionSettings(valley_min_canonical_repeats=5)
    assert _valley_pair(tmp_path, valleys, 24, settings=settings) == [6, 12]


def test_valley_ranking_matches_former_key_at_default_layout() -> None:
    valleys = [(6, 10.0), (12, 40.0), (20, 60.0), (41, 72.0), (69, 110.0)]
    fixed = DEFAULT_SETTINGS.reference_layout.fixed_repeat_count
    old = sorted(valleys, key=lambda x: x[1] / ((x[0] + 9) * 60))
    new = sorted(valleys, key=lambda x: x[1] / (x[0] + fixed))
    assert old == new


def test_valley_ranking_uses_reference_layout_fixed_count(tmp_path: Path) -> None:
    valleys = {10: 10, 40: 30, 70: 45}
    assert _valley_pair(tmp_path, valleys, 74) == [10, 70]
    ids = tuple(str(i) for i in range(200))
    layout = ReferenceLayoutSettings(pre=ids[:100], after=ids[100:])
    assert _valley_pair(tmp_path, valleys, 74, reference_layout=layout) == [10, 40]


# split_cluster_by_read_length ---------------------------------------------------------------


def _length_lines(lengths: dict[int, int]) -> list[str]:
    return [
        _sam(f"r{length}_{i}", "contig_31", cigar=f"{length}M", seq="A" * length)
        for length, count in lengths.items()
        for i in range(count)
    ]


CLUSTER = {"center": 32, "total_reads": 40, "contigs": [(30, 5), (31, 15), (32, 10), (33, 10)]}
BIMODAL = {2430: 20, 2550: 20}


def _centres(tmp_path: Path, lengths: dict[int, int], **kwargs: Any) -> list[int] | None:
    bam = tmp_path / "m.bam"
    bam.touch()
    lines = _length_lines(lengths)
    result = split_cluster_by_read_length(
        bam, CLUSTER, run_tool_iter_func=lambda cmd: lines, **kwargs
    )
    return None if result is None else [sub["center"] for sub in result]


def test_read_length_split_uses_reference_layout_fixed_count(tmp_path: Path) -> None:
    """The former literal 9 fixed repeats now follows the configured layout."""
    assert _centres(tmp_path, BIMODAL) == [31, 33]
    layout = ReferenceLayoutSettings(pre=("1", "2", "3", "4"), after=("6", "7", "8", "9"))
    assert _centres(tmp_path, BIMODAL, reference_layout=layout) == [32, 34]


def test_read_length_split_default_unit_is_the_dictionary_unit(tmp_path: Path) -> None:
    assert _centres(tmp_path, BIMODAL, repeat_length_bp=60) == _centres(tmp_path, BIMODAL)


def test_read_length_split_honours_a_non_sixty_unit(tmp_path: Path) -> None:
    lengths = {2030: 20, 2130: 20}  # 100 bp apart: two 50 bp units, not a whole 60 bp unit
    assert _centres(tmp_path, lengths) is None
    assert _centres(tmp_path, lengths, repeat_length_bp=50) == [31, 33]


@pytest.mark.parametrize(
    ("lengths", "values", "expected"),
    [
        pytest.param(BIMODAL, {"read_length_split_min_reads": 21}, None, id="min_reads"),
        pytest.param(BIMODAL, {"read_length_split_min_fraction": 0.9}, None, id="min_fraction"),
        pytest.param(BIMODAL, {"read_length_split_bin_bp": 200}, None, id="bin_bp"),
        pytest.param(BIMODAL, {"read_length_split_min_delta_bp": 121}, None, id="min_delta_bp"),
        pytest.param(BIMODAL, {"read_length_split_max_delta_bp": 119}, None, id="max_delta_bp"),
        pytest.param(BIMODAL, {"read_length_split_offset_bp": 90}, [30, 32], id="offset_bp"),
        pytest.param(
            {2430: 20, 2560: 20},
            {"read_length_split_unit_tolerance_bp": 5},
            None,
            id="unit_tolerance_bp",
        ),
    ],
)
def test_read_length_split_setting_is_honoured(
    tmp_path: Path, lengths: dict[int, int], values: dict[str, Any], expected: list[int] | None
) -> None:
    assert _centres(tmp_path, lengths) is not None
    settings = AlleleSelectionSettings(**values)
    assert _centres(tmp_path, lengths, settings=settings) == expected


def test_read_length_split_keyword_overrides_take_precedence(tmp_path: Path) -> None:
    assert _centres(tmp_path, BIMODAL, min_reads=21) is None
    assert _centres(tmp_path, BIMODAL, min_fraction=0.9) is None


def test_read_length_split_rejects_tolerance_of_half_a_unit(tmp_path: Path) -> None:
    settings = AlleleSelectionSettings(read_length_split_unit_tolerance_bp=30)
    with pytest.raises(ValueError, match="read_length_split_unit_tolerance_bp"):
        split_cluster_by_read_length(
            tmp_path / "missing.bam", CLUSTER, settings=settings, repeat_length_bp=60
        )


# _build_allele_info: ONT refinement shift -------------------------------------------------


@pytest.mark.parametrize(("min_shift_ont", "canonical"), [(2, 52), (0, 50)])
def test_ont_refinement_shift_floor_is_configurable(min_shift_ont: int, canonical: int) -> None:
    cluster = {"center": 50, "total_reads": 30, "contigs": [(50, 15), (52, 15)]}
    settings = AlleleSelectionSettings(refinement_min_shift_ont=min_shift_ont)
    info = _build_allele_info(cluster, "contig_52", settings=settings, platform="ont")
    assert info["canonical_repeats"] == canonical


# detect_alleles: minority floor and wiring -------------------------------------------------


def _detect(tmp_path: Path, settings: AlleleSelectionSettings | None = None) -> tuple[dict, Any]:
    metrics = {
        "contig_50": {"primary_alignment_records": 100},
        "contig_70": {"primary_alignment_records": 3},
    }

    def refine(bam: Path, contigs: list[str], **kwargs: Any) -> dict:
        return {"best_contig": contigs[0], "metrics": {c: metrics[c] for c in contigs}}

    valid = DominanceScore("contig_50", "contig_70", 0, 3, 0, 3, True)
    with (
        patch("muc_one_span.alleles.refine_peak_contig", side_effect=refine),
        patch("muc_one_span.alleles.split_cluster_by_read_length", return_value=None) as by_len,
        patch("muc_one_span.alleles._split_cluster_by_indel", return_value=None) as by_indel,
        patch("muc_one_span.alleles.extract_read_scores_for_contigs", return_value={"r": {}}),
        patch(
            "muc_one_span.alleles.evaluate_candidate_pair_dominance", return_value=valid
        ) as dominance,
    ):
        result = detect_alleles(
            {50: 100, 70: 3},
            10,
            bam_path=tmp_path / "missing.bam",
            settings=settings,
            repeat_length_bp=55,
        )
    return result, (by_len, by_indel, dominance)


def test_minority_alignment_floor_is_configurable(tmp_path: Path) -> None:
    result, _ = _detect(tmp_path)
    assert result["allele_2"]["contig_name"] == "contig_70"
    settings = AlleleSelectionSettings(minority_min_alignment_records=4)
    result, _ = _detect(tmp_path, settings)
    assert result["same_length"] is True


def test_detect_alleles_forwards_settings_layout_and_unit(tmp_path: Path) -> None:
    settings = AlleleSelectionSettings(
        dominance_close_candidate_repeats=4, dominance_zero_primary_extra_reads=3
    )
    _, (by_len, by_indel, dominance) = _detect(tmp_path, settings)
    layout = DEFAULT_SETTINGS.reference_layout
    assert by_len.call_args.kwargs["settings"] is settings
    assert by_len.call_args.kwargs["reference_layout"] == layout
    assert by_len.call_args.kwargs["repeat_length_bp"] == 55
    assert by_indel.call_args.kwargs["settings"] is settings
    assert by_indel.call_args.kwargs["reference_layout"] == layout
    assert dominance.call_args.kwargs["close_candidate_repeats"] == 4
    assert dominance.call_args.kwargs["zero_primary_extra_reads"] == 3


# evaluate_candidate_pair_dominance ---------------------------------------------------------


def test_dominance_defaults_come_from_settings() -> None:
    assert DEFAULT_DELTA_AS["hifi"] == DEFAULT.score_margin_hifi
    assert DEFAULT_DELTA_AS["ont"] == DEFAULT.score_margin_ont
    assert DEFAULT_MIN_DOMINANT_READS["hifi"] == DEFAULT.min_dominant_reads_hifi
    assert DEFAULT_MIN_DOMINANT_READS["ont"] == DEFAULT.min_dominant_reads_ont
    assert DEFAULT.min_dominance_ratio == DEFAULT_MIN_RATIO


@pytest.mark.parametrize(("close", "valid"), [(None, False), (5, True)])
def test_close_candidate_repeats_is_configurable(close: int | None, valid: bool) -> None:
    scores = {f"r{i}": {"contig_55": 1000} for i in range(10)}
    result = evaluate_candidate_pair_dominance(
        scores, "contig_50", "contig_55", "hifi", close_candidate_repeats=close
    )
    assert result.is_valid_second_allele is valid


@pytest.mark.parametrize(("extra", "valid"), [(None, False), (1, True)])
def test_zero_primary_extra_reads_is_configurable(extra: int | None, valid: bool) -> None:
    d2 = DEFAULT.min_dominant_reads_hifi + 1
    scores = {f"r{i}": {"contig_50": 900, "contig_60": 1000} for i in range(d2)}
    result = evaluate_candidate_pair_dominance(
        scores,
        "contig_50",
        "contig_60",
        "hifi",
        c2_primary_records=0,
        zero_primary_extra_reads=extra,
    )
    assert result.is_valid_second_allele is valid
