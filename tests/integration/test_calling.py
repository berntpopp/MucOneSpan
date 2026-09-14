"""Integration tests for allele extraction, VCF filtering, and Clair3."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from muc_one_span.calling import extract_allele_reads, run_clair3
from muc_one_span.vcf import filter_vcf
from tests.conftest import requires_bcftools, requires_clair3, requires_samtools


@requires_samtools
@pytest.mark.integration
class TestExtractAlleleReads:
    @pytest.mark.parametrize(
        "contigs,expected",
        [("contig_1", {"read1", "read2"}), (["contig_1", "contig_2"], {"read1", "read2", "read3"})],
    )
    def test_extracts_only_requested_contigs(
        self, tmp_path: Path, aligned_bam: Path, contigs: str | list[str], expected: set[str]
    ) -> None:
        bam = extract_allele_reads(aligned_bam, contigs, tmp_path / "extracted")
        assert Path(str(bam) + ".bai").is_file()
        subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
        records = subprocess.check_output(["samtools", "view", str(bam)], text=True)
        assert {line.split("\t")[0] for line in records.splitlines()} == expected
        assert len(records.splitlines()) == len(expected)
        # Query via the index, which must also be usable after subsetting.
        indexed = subprocess.check_output(
            ["samtools", "view", "-c", str(bam), "contig_1"], text=True
        )
        assert int(indexed) == 2


@requires_bcftools
@requires_samtools
@pytest.mark.integration
class TestFilterVcf:
    def test_filters_qual_and_pass_without_info_depth(
        self,
        tmp_path: Path,
        reference_data: tuple[Path, dict[str, str]],
        variant_data: tuple[Path, dict[str, str]],
    ) -> None:
        reference, _ = reference_data
        vcf, _ = variant_data
        filtered = filter_vcf(vcf, reference, tmp_path / "filtered", min_qual=5, min_dp=100)
        assert Path(str(filtered) + ".csi").is_file()
        query = subprocess.check_output(
            ["bcftools", "query", "-r", "contig_1", "-f", "%POS\t%QUAL\t%FILTER\n", str(filtered)],
            text=True,
        )
        assert query.splitlines() == ["500\t30\tPASS", "1000\t30\tPASS"]
        assert not (filtered.parent / "normalized.vcf.gz").exists()

    def test_header_only_vcf_remains_valid(
        self,
        tmp_path: Path,
        reference_data: tuple[Path, dict[str, str]],
        variant_data: tuple[Path, dict[str, str]],
    ) -> None:
        reference, _ = reference_data
        vcf, _ = variant_data
        vcf.write_text(
            "".join(line + "\n" for line in vcf.read_text().splitlines() if line.startswith("#"))
        )
        filtered = filter_vcf(vcf, reference, tmp_path / "empty", min_qual=5)
        records = subprocess.check_output(["bcftools", "view", "-H", str(filtered)], text=True)
        assert records == ""
        assert Path(str(filtered) + ".csi").is_file()


@requires_clair3
@requires_samtools
@requires_bcftools
@pytest.mark.integration
class TestRunClair3:
    def test_calls_supported_snp(
        self, tmp_path: Path, reference_data: tuple[Path, dict[str, str]], clair3_model: Path
    ) -> None:
        """Forty perfect HiFi reads support a single alternate base."""
        reference, sequences = reference_data
        seq = sequences["contig_1"]
        alternate = next(base for base in "ACGT" if base != seq[1199])
        read = seq[:1199] + alternate + seq[1200:]
        reference.write_text(f">contig_1\n{seq}\n")
        subprocess.run(["samtools", "faidx", str(reference)], check=True)
        sam = tmp_path / "clair3.sam"
        sam.write_text(
            f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:contig_1\tLN:{len(seq)}\n"
            + "".join(
                f"read{i}\t0\tcontig_1\t1\t60\t{len(seq)}M\t*\t0\t0\t{read}\t{'I' * len(read)}\n"
                for i in range(40)
            )
        )
        bam = tmp_path / "clair3.bam"
        subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
        subprocess.run(["samtools", "index", str(bam)], check=True)
        vcf = run_clair3(
            bam, reference, tmp_path / "calls", model_path=str(clair3_model), threads=2
        )
        records = subprocess.check_output(
            ["bcftools", "query", "-f", "%POS\t%REF\t%ALT\n", str(vcf)], text=True
        )
        assert f"1200\t{seq[1199]}\t{alternate}" in records.splitlines()
