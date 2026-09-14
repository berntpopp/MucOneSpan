"""Distinguish reference fit from record support and avoid unsupported splits."""

from pathlib import Path
from unittest.mock import patch

from muc_one_span.alleles import _build_allele_info, _split_cluster_by_indel, refine_peak_contig


def sam(name: str, contig: int, flag: int = 0, cigar: str = "60M") -> str:
    return f"{name}\t{flag}\tcontig_{contig}\t1\t60\t{cigar}\t*\t0\t0\tAAAA\tIIII\tAS:i:60"


def test_flat_indel_plateau_is_not_two_valleys(tmp_path: Path) -> None:
    cluster = {"center": 43, "total_reads": 30, "contigs": [(40, 10), (43, 10), (46, 10)]}
    with patch(
        "muc_one_span.alleles.run_tool_iter", return_value=iter(sam("r", c) for c in (40, 43, 46))
    ):
        assert _split_cluster_by_indel(tmp_path / "mapping.bam", cluster) is None


def test_primary_records_preserve_qname_collisions_and_secondary_fit(tmp_path: Path) -> None:
    lines = [
        sam("collision", 51),
        sam("collision", 51),
        sam("collision", 51, 256),
        sam("r", 51, 2048),
    ]
    with patch("muc_one_span.alleles.run_tool_iter", return_value=iter(lines)):
        result = refine_peak_contig(tmp_path / "mapping.bam", ["contig_51"])
    metric = result["metrics"]["contig_51"]
    assert metric["reads"] == 4
    assert metric["primary_alignment_records"] == 2
    assert metric["secondary_alignment_records"] == 1
    assert metric["supplementary_alignment_records"] == 1
    assert result["best_contig"] == "contig_51"


def test_candidate_count_and_selected_reference_length_are_explicit() -> None:
    result = _build_allele_info(
        {"center": 72, "total_reads": 27, "contigs": [(72, 27)]}, "contig_69"
    )
    assert result["length"] == 81
    assert result["reference_length"] == 78
    assert result["alignment_records"] == 27
    assert result["molecule_count"] is None
    assert result["support_basis"] == "alignment_records"


def test_detected_support_is_from_primary_records_not_secondary_fit(tmp_path: Path) -> None:
    from muc_one_span.alleles import detect_alleles

    lines = [sam("collision", 51), sam("collision", 51), sam("collision", 51, 256)]
    with patch("muc_one_span.alleles.run_tool_iter", side_effect=lambda _: iter(lines)):
        result = detect_alleles({51: 30}, bam_path=tmp_path / "mapping.bam")
    assert result["allele_1"]["primary_alignment_records"] == 2
    assert result["allele_1"]["alignment_records"] == 30
    assert result["observed_length_candidates"] == 1
    assert result["allele_2"]["candidate_duplicate_of"] == "allele_1"


def test_length_selection_evidence_retains_dropped_and_unselected_candidates() -> None:
    from muc_one_span.alleles import detect_alleles

    result = detect_alleles({50: 9, 51: 10, 80: 11, 110: 10}, min_coverage=10)

    assert sorted((result["allele_1"]["length"], result["allele_2"]["length"])) == [60, 89]
    for allele in (result["allele_1"], result["allele_2"]):
        evidence = allele["length_selection_evidence"]
        assert evidence["minimum_coverage"] == 10
        assert evidence["support_unit"] == "alignment_records_not_molecules"
        assert evidence["molecule_count"] is None
        assert evidence["excluded_subthreshold_alignment_records"] == 9
        assert evidence["excluded_subthreshold_primary_alignment_records"] is None
        assert evidence["excluded_subthreshold_contigs"] == [
            {
                "contig_name": "contig_50",
                "alignment_records": 9,
                "primary_alignment_records": None,
                "molecule_count": None,
            }
        ]
        assert evidence["unselected_passing_clusters"] == [
            {
                "center": 110,
                "alignment_records": 10,
                "contigs": [{"contig_name": "contig_110", "alignment_records": 10}],
            }
        ]


def test_excluded_primary_records_count_colliding_qnames_in_one_query(tmp_path: Path) -> None:
    from muc_one_span.alleles import detect_alleles

    bam = tmp_path / "mapping.bam"
    bam.touch()
    commands: list[list[str]] = []

    def records(args: list[str]):
        commands.append(args)
        contigs = args[3:]
        if contigs == ["contig_100", "contig_101"]:
            return iter(
                [
                    sam("collision", 100),
                    sam("collision", 100),
                    sam("collision", 100, 256),
                    sam("other", 101),
                    sam("unmapped", 101, 4),
                ]
            )
        return iter(sam("selected", int(contig.removeprefix("contig_"))) for contig in contigs)

    with patch("muc_one_span.alleles.run_tool_iter", side_effect=records):
        result = detect_alleles({51: 30, 80: 20, 100: 4, 101: 3}, min_coverage=10, bam_path=bam)

    evidence = result["allele_1"]["length_selection_evidence"]
    assert evidence["excluded_subthreshold_primary_alignment_records"] == 3
    assert [
        row["primary_alignment_records"] for row in evidence["excluded_subthreshold_contigs"]
    ] == [
        2,
        1,
    ]
    assert sum(command[3:] == ["contig_100", "contig_101"] for command in commands) == 1


def test_empty_excluded_evidence_has_zero_counts_without_bam() -> None:
    from muc_one_span.alleles import detect_alleles

    result = detect_alleles({51: 30}, min_coverage=10)
    evidence = result["allele_1"]["length_selection_evidence"]
    assert evidence["excluded_subthreshold_contigs"] == []
    assert evidence["excluded_subthreshold_alignment_records"] == 0
    assert evidence["excluded_subthreshold_primary_alignment_records"] == 0


def test_gap_setting_changes_observed_length_candidates() -> None:
    from muc_one_span.alleles import detect_alleles
    from muc_one_span.settings import AlleleSelectionSettings

    counts = {40: 20, 44: 20}
    assert detect_alleles(counts)["observed_length_candidates"] == 1
    result = detect_alleles(counts, settings=AlleleSelectionSettings(min_gap=4))
    assert result["observed_length_candidates"] == 2
    assert [result[key]["length"] for key in ("allele_1", "allele_2")] == [49, 53]


def test_refinement_and_fixed_layout_count_are_configurable() -> None:
    from muc_one_span.alleles import detect_alleles
    from muc_one_span.settings import AlleleSelectionSettings, ReferenceLayoutSettings

    cluster = {"center": 40, "total_reads": 20, "contigs": [(40, 20)]}
    layout = ReferenceLayoutSettings(pre=("1", "2"), after=("9",))
    result = _build_allele_info(
        cluster,
        "contig_42",
        reference_layout=layout,
        settings=AlleleSelectionSettings(refinement_max_shift=2),
    )
    assert result["length"] == result["reference_length"] == 45
    assert result["fixed_repeat_count"] == 3
    default = _build_allele_info(cluster, "contig_42")
    assert default["length"] == 49
    assert default["fixed_repeat_count"] == 9
    aliases = detect_alleles({40: 20}, reference_layout=layout)
    for key in ("allele_1", "allele_2"):
        assert aliases[key]["length"] == 43
        assert aliases[key]["fixed_repeat_count"] == 3


def test_valley_settings_control_points_and_separation(tmp_path: Path) -> None:
    from muc_one_span.settings import AlleleSelectionSettings

    cluster = {"center": 41, "total_reads": 30, "contigs": [(40, 10), (41, 10), (42, 10)]}
    lines = [sam("r", 40), sam("r", 41, cigar="50M10I"), sam("r", 42)]
    with patch("muc_one_span.alleles.run_tool_iter", side_effect=lambda _: iter(lines)):
        assert _split_cluster_by_indel(tmp_path / "b.bam", cluster) is None
        result = _split_cluster_by_indel(
            tmp_path / "b.bam", cluster, settings=AlleleSelectionSettings(valley_min_separation=2)
        )
        assert result is not None
        assert [item["center"] for item in result] == [40, 42]
        assert (
            _split_cluster_by_indel(
                tmp_path / "b.bam",
                cluster,
                settings=AlleleSelectionSettings(valley_min_points=4, valley_min_separation=2),
            )
            is None
        )
