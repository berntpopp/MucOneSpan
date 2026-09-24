"""Tests for genotype and phase handling in distinct-length calling."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.calling import call_variants_per_allele
from muc_one_span.settings import CallingSettings


def _make_alleles():
    return {
        "homozygous": False,
        "same_length": False,
        "allele_1": {
            "length": 37,
            "reads": 1000,
            "canonical_repeats": 28,
            "contig_name": "contig_28",
            "cluster_contigs": ["contig_28"],
        },
        "allele_2": {
            "length": 70,
            "reads": 800,
            "canonical_repeats": 61,
            "contig_name": "contig_61",
            "cluster_contigs": ["contig_61"],
        },
    }


@patch("muc_one_span.calling._extract_and_remap_reads")
@patch("muc_one_span.calling.run_clair3")
@patch("muc_one_span.calling.filter_vcf")
def test_distinct_length_unphased_retains_iupac_selector(
    mock_filter, mock_clair3, mock_remap, tmp_path
):
    """Unphased multi-het candidate must receive 'I' selector without fabricated phase credit."""
    mock_remap.return_value = tmp_path / "mapped.bam"
    mock_clair3.return_value = tmp_path / "raw.vcf.gz"
    vcf1 = tmp_path / "allele_1.vcf.gz"
    vcf2 = tmp_path / "allele_2.vcf.gz"
    mock_filter.side_effect = [vcf1, vcf2]

    unphased_variants = [
        {
            "chrom": "contig_28",
            "pos": 100,
            "ref": "G",
            "alt": "GC",
            "genotype": "0/1",
            "sample": "S1",
        },
        {
            "chrom": "contig_28",
            "pos": 200,
            "ref": "C",
            "alt": "A",
            "genotype": "0/1",
            "sample": "S1",
        },
    ]
    homozygous_variants = [
        {
            "chrom": "contig_61",
            "pos": 150,
            "ref": "A",
            "alt": "G",
            "genotype": "1/1",
            "sample": "S1",
        },
    ]

    alleles = _make_alleles()
    with patch(
        "muc_one_span.calling.parse_vcf_genotypes",
        side_effect=[unphased_variants, homozygous_variants],
    ):
        result = call_variants_per_allele(
            tmp_path / "in.bam",
            tmp_path / "ref.fa",
            alleles,
            tmp_path / "out",
        )

    assert result == {"allele_1": vcf1, "allele_2": vcf2}
    a1 = alleles["allele_1"]
    assert a1["phase_status"] == "unphased"
    assert a1["consensus_haplotype"] == "I"
    assert a1["consensus_policy"] == "genotype_iupac_candidate"
    assert a1["sequence_source"].endswith(":GTI")
    assert a1["independent_haplotype_evidence"] is False

    a2 = alleles["allele_2"]
    assert a2["consensus_haplotype"] == 1
    assert a2["independent_haplotype_evidence"] is True


def _variant(
    pos: int, ref: str, alt: str, genotype: str, phase_set: str | None = None
) -> dict[str, Any]:
    return {
        "chrom": "contig_71",
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "qual": 29.93,
        "genotype": genotype,
        "phase_set": phase_set,
        "sample": "S1",
    }


def _run_distinct(
    tmp_path: Path, allele_2_variants: list[dict], **call_kwargs: Any
) -> tuple[dict, list[dict[str, Any]]]:
    """Run distinct-length calling with tools mocked; allele 2 is MP4-shaped (80 units)."""
    alleles = _make_alleles()
    alleles["allele_2"].update(
        length=80, canonical_repeats=71, contig_name="contig_71", cluster_contigs=["contig_71"]
    )
    by_allele: dict[str, list[dict]] = {"allele_1": [], "allele_2": allele_2_variants}
    filter_calls: list[dict[str, Any]] = []

    def fake_filter(vcf: Path, reference: Path, out_dir: Path, **kwargs: Any) -> Path:
        filter_calls.append(kwargs)
        return out_dir / "variants.vcf.gz"

    with (
        patch("muc_one_span.calling._extract_and_remap_reads", return_value=tmp_path / "m.bam"),
        patch("muc_one_span.calling.run_clair3", return_value=tmp_path / "raw.vcf.gz"),
        patch("muc_one_span.calling.filter_vcf", side_effect=fake_filter),
        patch(
            "muc_one_span.calling.parse_vcf_genotypes",
            side_effect=lambda path, sample=None: by_allele[path.parent.name],
        ),
    ):
        call_variants_per_allele(
            tmp_path / "in.bam", tmp_path / "ref.fa", alleles, tmp_path / "out", **call_kwargs
        )
    return alleles, filter_calls


def test_mp4_shaped_single_het_insertion_is_not_resolved_as_reference(tmp_path: Path) -> None:
    """ERR15277569 shape: a 0/1 G>GC left after the allele-fraction rule must not become REF."""
    alleles, _ = _run_distinct(tmp_path, [_variant(3432, "G", "GC", "0/1")])
    a2 = alleles["allele_2"]
    assert a2["phase_status"] == "single_heterozygous_unordered"
    assert a2["consensus_haplotype"] == "I"
    assert a2["consensus_policy"] == "genotype_iupac_candidate"
    assert a2["allele_genotype_status"] == "heterozygous_within_length_partition"
    assert a2["heterozygous_sites"] == [
        {"chrom": "contig_71", "pos": 3432, "ref": "G", "alt": "GC", "genotype": "0/1"}
    ]
    assert a2["independent_haplotype_evidence"] is False
    assert alleles["allele_1"]["consensus_haplotype"] == 1
    assert alleles["allele_1"]["independent_haplotype_evidence"] is True


@pytest.mark.parametrize(
    ("variants", "selector", "status", "independent"),
    [
        pytest.param(
            [_variant(3432, "G", "GC", "1/1")],
            1,
            "allele_specific_resolved",
            True,
            id="homozygous_alt",
        ),
        pytest.param([], 1, "allele_specific_resolved", True, id="no_retained_variants"),
        pytest.param(
            [_variant(3432, "G", "GC", "0/1")],
            "I",
            "heterozygous_within_length_partition",
            False,
            id="single_site_het",
        ),
        pytest.param(
            [_variant(100, "G", "GC", "0|1", "100"), _variant(200, "C", "A", "1|0", "100")],
            "I",
            "heterozygous_within_length_partition",
            False,
            id="phased",
        ),
        pytest.param(
            [_variant(100, "G", "GC", "0/1"), _variant(200, "C", "A", "0/1")],
            "I",
            "heterozygous_within_length_partition",
            False,
            id="multi_site_unphased",
        ),
        pytest.param(
            [_variant(100, "GC", "G", "1/1"), _variant(101, "C", "A", "1/1")],
            "I",
            "unresolved_genotype_records",
            False,
            id="conflicting_records",
        ),
    ],
)
def test_length_partition_selector_follows_allele_specific_genotypes(
    tmp_path: Path, variants: list[dict], selector: int | str, status: str, independent: bool
) -> None:
    alleles, _ = _run_distinct(tmp_path, variants)
    a2 = alleles["allele_2"]
    assert a2["consensus_haplotype"] == selector
    assert a2["allele_genotype_status"] == status
    assert a2["independent_haplotype_evidence"] is independent


def test_distinct_path_forwards_allele_fraction_settings(tmp_path: Path) -> None:
    settings = CallingSettings(haploid_alt_fraction=0.6, haploid_ref_fraction=0.1)
    _, filter_calls = _run_distinct(tmp_path, [], settings=settings)
    assert len(filter_calls) == 2
    for kwargs in filter_calls:
        assert (kwargs["haploid_alt_fraction"], kwargs["haploid_ref_fraction"]) == (0.6, 0.1)


def test_default_haploid_filter_is_recorded(tmp_path: Path) -> None:
    alleles, filter_calls = _run_distinct(tmp_path, [])
    assert {call["haploid_min_qual"] for call in filter_calls} == {4.0}
    assert {call["haploid_majority"] for call in filter_calls} == {True}
    assert alleles["allele_2"]["variant_filter"] == {
        "min_qual": 4.0,
        "min_qual_source": "calling.haploid_min_qual",
        "haploid_majority": True,
    }


def test_null_haploid_min_qual_follows_min_qual(tmp_path: Path) -> None:
    settings = CallingSettings(haploid_min_qual=None, haploid_majority=False)
    alleles, filter_calls = _run_distinct(tmp_path, [], min_qual=12.0, settings=settings)
    assert {call["haploid_min_qual"] for call in filter_calls} == {12.0}
    assert {call["haploid_majority"] for call in filter_calls} == {False}
    assert alleles["allele_1"]["variant_filter"] == {
        "min_qual": 12.0,
        "min_qual_source": "run.min_qual",
        "haploid_majority": False,
    }


def test_overridden_explicit_min_qual_is_logged(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    _, filter_calls = _run_distinct(tmp_path, [], min_qual=12.0)
    assert {call["haploid_min_qual"] for call in filter_calls} == {4.0}
    assert "not applied to length-partitioned calls" in caplog.text
