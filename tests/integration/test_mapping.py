"""Integration tests for read mapping and BAM conversion with real tools."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from open_pacmuci.mapping import bam_to_fastq, map_reads
from tests.conftest import requires_minimap2, requires_samtools


@requires_minimap2
@requires_samtools
@pytest.mark.integration
class TestMapReads:
    @pytest.mark.parametrize("input_format", ["fastq", "bam"])
    def test_maps_reads_to_sorted_indexed_bam(
        self,
        tmp_path: Path,
        reference_data: tuple[Path, dict[str, str]],
        aligned_bam: Path,
        input_format: str,
    ) -> None:
        reference, sequences = reference_data
        if input_format == "bam":
            reads = aligned_bam
        else:
            reads = tmp_path / "reads.fq"
            reads.write_text(
                "".join(
                    f"@{name}\n{sequences[contig][start : start + 1500]}\n+\n{'I' * 1500}\n"
                    for name, contig, start in (
                        ("read1", "contig_1", 100),
                        ("read2", "contig_1", 500),
                        ("read3", "contig_2", 200),
                    )
                )
            )
        bam = map_reads(reads, reference, tmp_path / "mapped", threads=1)
        assert bam.name == "mapping.bam"
        assert Path(str(bam) + ".bai").is_file()
        subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
        header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
        assert "SO:coordinate" in header
        records = subprocess.check_output(["samtools", "view", "-F", "2308", str(bam)], text=True)
        fields = [line.split("\t") for line in records.splitlines()]
        assert [(r[0], r[2], int(r[3]), r[5]) for r in fields] == [
            ("read1", "contig_1", 101, "1500M"),
            ("read2", "contig_1", 501, "1500M"),
            ("read3", "contig_2", 201, "1500M"),
        ]
        assert all(int(r[4]) >= 50 for r in fields)


@requires_samtools
@pytest.mark.integration
class TestBamToFastq:
    def test_preserves_read_names_sequences_and_qualities(
        self, tmp_path: Path, reference_data: tuple[Path, dict[str, str]], aligned_bam: Path
    ) -> None:
        _, sequences = reference_data
        fastq = bam_to_fastq(aligned_bam, tmp_path / "converted")
        lines = fastq.read_text().splitlines()
        assert len(lines) == 12
        records = {lines[i]: (lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)}
        assert records == {
            "@read1": (sequences["contig_1"][100:1600], "I" * 1500),
            "@read2": (sequences["contig_1"][500:2000], "I" * 1500),
            "@read3": (sequences["contig_2"][200:1700], "I" * 1500),
        }
