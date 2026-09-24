"""Tests for length-partition selection on the read-phased (haplotagged) same-length path."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.calling import disambiguate_same_length_alleles
from muc_one_span.settings import CallingSettings


def _variant(pos: int, ref: str, alt: str, genotype: str) -> dict[str, Any]:
    return {
        "chrom": "contig_60",
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "qual": 29.93,
        "genotype": genotype,
        "phase_set": None,
        "sample": "S1",
    }


def _run_haplotagged(
    tmp_path: Path, hp1: list[dict], hp2: list[dict], **call_kwargs: Any
) -> tuple[dict, list[dict[str, Any]]]:
    """Run same-length calling through the haplotag split with all tools mocked."""
    alleles = {
        "same_length": True,
        "allele_1": {"contig_name": "contig_60", "cluster_contigs": ["contig_60"], "length": 60},
    }
    by_dir = {"merged": [], "allele_1": hp1, "allele_2": hp2}
    filter_calls: list[dict[str, Any]] = []

    def fake_filter(vcf: Path, reference: Path, out_dir: Path, **kwargs: Any) -> Path:
        filter_calls.append({"dir": out_dir.name, **kwargs})
        return out_dir / "variants.vcf.gz"

    evidence = {"status": "phased", "output_phase_status": "phased"}
    with (
        patch("muc_one_span.calling._extract_and_remap_reads", return_value=tmp_path / "m.bam"),
        patch("muc_one_span.calling.run_clair3", return_value=tmp_path / "raw.vcf.gz"),
        patch("muc_one_span.calling.filter_vcf", side_effect=fake_filter),
        patch(
            "muc_one_span.calling.phase_same_length_reads",
            side_effect=lambda vcf, *args, **kwargs: (vcf, evidence),
        ),
        patch(
            "muc_one_span.calling.haplotag_and_split_reads",
            return_value=(tmp_path / "hp1.bam", tmp_path / "hp2.bam", 25, 20),
        ),
        patch(
            "muc_one_span.calling.parse_vcf_genotypes",
            side_effect=lambda path, sample=None: by_dir[path.parent.name],
        ),
    ):
        disambiguate_same_length_alleles(
            tmp_path / "in.bam",
            tmp_path / "ref.fa",
            alleles,
            tmp_path / "out",
            read_phase=True,
            **call_kwargs,
        )
    return alleles, filter_calls


def test_haplotag_het_record_in_ambiguous_band_is_not_resolved_as_reference(
    tmp_path: Path,
) -> None:
    """A 0/1 left after the AD-fraction rule must select 'I', never GT1 (= REF)."""
    alleles, _ = _run_haplotagged(
        tmp_path, [_variant(3432, "G", "GC", "0/1")], [_variant(150, "A", "G", "1/1")]
    )
    a1 = alleles["allele_1"]
    assert a1["consensus_haplotype"] == "I"
    assert a1["consensus_policy"] == "genotype_iupac_candidate"
    assert a1["allele_genotype_status"] == "heterozygous_within_length_partition"
    assert a1["heterozygous_sites"] == [
        {"chrom": "contig_60", "pos": 3432, "ref": "G", "alt": "GC", "genotype": "0/1"}
    ]
    assert a1["independent_haplotype_evidence"] is False
    a2 = alleles["allele_2"]
    assert a2["consensus_haplotype"] == 1
    assert a2["allele_genotype_status"] == "allele_specific_resolved"
    assert a2["independent_haplotype_evidence"] is True


@pytest.mark.parametrize(
    ("variants", "selector", "status", "independent"),
    [
        pytest.param(
            [_variant(3432, "G", "GC", "1/1")], 1, "allele_specific_resolved", True, id="hom_alt"
        ),
        pytest.param([], 1, "allele_specific_resolved", True, id="no_retained_variants"),
        pytest.param(
            [_variant(100, "G", "GC", "0/1"), _variant(200, "C", "A", "0/1")],
            "I",
            "heterozygous_within_length_partition",
            False,
            id="multi_site_het",
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
def test_haplotag_selector_follows_allele_specific_genotypes(
    tmp_path: Path, variants: list[dict], selector: int | str, status: str, independent: bool
) -> None:
    alleles, _ = _run_haplotagged(tmp_path, variants, [])
    a1 = alleles["allele_1"]
    assert a1["consensus_haplotype"] == selector
    assert a1["allele_genotype_status"] == status
    assert a1["independent_haplotype_evidence"] is independent


def test_haplotag_path_forwards_filter_settings_and_records_provenance(tmp_path: Path) -> None:
    settings = CallingSettings(
        haploid_alt_fraction=0.6, haploid_ref_fraction=0.1, haploid_majority=False
    )
    alleles, filter_calls = _run_haplotagged(tmp_path, [], [], settings=settings)
    haplotype_calls = [call for call in filter_calls if call["dir"] != "merged"]
    assert len(haplotype_calls) == 2
    for kwargs in haplotype_calls:
        assert (kwargs["haploid_alt_fraction"], kwargs["haploid_ref_fraction"]) == (0.6, 0.1)
        assert kwargs["haploid_majority"] is False
        assert kwargs["haploid_min_qual"] == 4.0
    expected = {
        "min_qual": 4.0,
        "min_qual_source": "calling.haploid_min_qual",
        "haploid_majority": False,
    }
    assert alleles["allele_1"]["variant_filter"] == expected
    assert alleles["allele_2"]["variant_filter"] == expected


def test_haplotag_diploid_genotypes_kept_by_configuration_are_unresolved(tmp_path: Path) -> None:
    """haploid_majority=false keeps 0/1 records; the haplotype must not drop the ALT."""
    alleles, _ = _run_haplotagged(
        tmp_path,
        [_variant(3432, "G", "GC", "0/1")],
        [_variant(3432, "G", "GC", "0/1")],
        settings=CallingSettings(haploid_majority=False),
    )
    for key in ("allele_1", "allele_2"):
        assert alleles[key]["consensus_haplotype"] == "I"
        assert alleles[key]["independent_haplotype_evidence"] is False
