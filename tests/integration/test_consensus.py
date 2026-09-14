"""Integration tests applying real bcftools variants to FASTA sequences."""

from __future__ import annotations

from pathlib import Path

import pytest

from open_pacmuci.consensus import build_consensus, build_consensus_per_allele
from open_pacmuci.vcf import filter_vcf
from tests.conftest import requires_bcftools, requires_samtools


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            sequences[name] = ""
        else:
            sequences[name] += line
    return sequences


@requires_bcftools
@requires_samtools
@pytest.mark.integration
class TestBuildConsensus:
    def test_applies_snp_and_insertion(
        self,
        tmp_path: Path,
        reference_data: tuple[Path, dict[str, str]],
        variant_data: tuple[Path, dict[str, str]],
    ) -> None:
        reference, _ = reference_data
        vcf, expected = variant_data
        filtered = filter_vcf(vcf, reference, tmp_path / "filtered", min_qual=5)
        consensus = build_consensus(reference, filtered, tmp_path / "consensus" / "result.fa")
        assert read_fasta(consensus) == expected


@requires_bcftools
@requires_samtools
@pytest.mark.integration
class TestBuildConsensusPerAllele:
    def test_extracts_correct_contigs_and_trims_flanks(
        self,
        tmp_path: Path,
        reference_data: tuple[Path, dict[str, str]],
        variant_data: tuple[Path, dict[str, str]],
    ) -> None:
        reference, _ = reference_data
        vcf, expected = variant_data
        filtered = filter_vcf(vcf, reference, tmp_path / "filtered", min_qual=5)
        alleles = {
            "allele_1": {"length": 10, "contig_name": "contig_1"},
            "allele_2": {"length": 11, "contig_name": "contig_2"},
        }
        outputs = build_consensus_per_allele(
            reference,
            dict.fromkeys(alleles, filtered),
            alleles,
            tmp_path / "consensus",
            flank_length=100,
        )
        assert set(outputs) == set(alleles)
        for allele, contig in (("allele_1", "contig_1"), ("allele_2", "contig_2")):
            assert read_fasta(outputs[allele]) == {f"{contig}_vntr": expected[contig][100:-100]}
