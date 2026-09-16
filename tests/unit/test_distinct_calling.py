"""Tests for genotype and phase handling in distinct-length calling."""

from __future__ import annotations

from unittest.mock import patch

from muc_one_span.calling import call_variants_per_allele


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


@patch("muc_one_span.calling._extract_and_remap_reads")
@patch("muc_one_span.calling.run_clair3")
@patch("muc_one_span.calling.filter_vcf")
def test_distinct_length_phased_and_single_het(mock_filter, mock_clair3, mock_remap, tmp_path):
    """Phased candidate and single-het candidate receive selector 1 and independent credit."""
    mock_remap.return_value = tmp_path / "mapped.bam"
    mock_clair3.return_value = tmp_path / "raw.vcf.gz"
    vcf1 = tmp_path / "allele_1.vcf.gz"
    vcf2 = tmp_path / "allele_2.vcf.gz"
    mock_filter.side_effect = [vcf1, vcf2]

    phased_variants = [
        {
            "chrom": "contig_28",
            "pos": 100,
            "ref": "G",
            "alt": "GC",
            "genotype": "0|1",
            "phase_set": "100",
            "sample": "S1",
        },
        {
            "chrom": "contig_28",
            "pos": 200,
            "ref": "C",
            "alt": "A",
            "genotype": "1|0",
            "phase_set": "100",
            "sample": "S1",
        },
    ]
    single_het = [
        {
            "chrom": "contig_61",
            "pos": 150,
            "ref": "A",
            "alt": "G",
            "genotype": "0/1",
            "sample": "S1",
        },
    ]

    alleles = _make_alleles()
    with patch(
        "muc_one_span.calling.parse_vcf_genotypes",
        side_effect=[phased_variants, single_het],
    ):
        result = call_variants_per_allele(
            tmp_path / "in.bam",
            tmp_path / "ref.fa",
            alleles,
            tmp_path / "out",
        )

    assert result == {"allele_1": vcf1, "allele_2": vcf2}
    a1 = alleles["allele_1"]
    assert a1["phase_status"] == "phased"
    assert a1["consensus_haplotype"] == 1
    assert a1["independent_haplotype_evidence"] is True

    a2 = alleles["allele_2"]
    assert a2["phase_status"] == "single_heterozygous_unordered"
    assert a2["consensus_haplotype"] == 1
    assert a2["independent_haplotype_evidence"] is True
