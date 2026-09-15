"""Mocked unit tests for mapping module (minimap2/samtools wrappers)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.mapping import bam_to_fastq, get_idxstats, map_reads


class TestBamToFastq:
    """Tests for bam_to_fastq."""

    @patch("muc_one_span.mapping.run_tool")
    def test_calls_samtools_fastq(self, mock_run_tool, tmp_path):
        """bam_to_fastq calls samtools fastq with the BAM path."""
        mock_run_tool.return_value = "@read1\nACGT\n+\nIIII\n"
        bam = tmp_path / "input.bam"
        bam.touch()

        bam_to_fastq(bam, tmp_path / "out")

        mock_run_tool.assert_called_once_with(["samtools", "fastq", str(bam)])

    @patch("muc_one_span.mapping.run_tool")
    def test_writes_fastq_output(self, mock_run_tool, tmp_path):
        """bam_to_fastq writes samtools stdout to extracted_reads.fq."""
        fastq_content = "@read1\nACGT\n+\nIIII\n"
        mock_run_tool.return_value = fastq_content
        bam = tmp_path / "input.bam"
        bam.touch()
        out_dir = tmp_path / "out"

        result = bam_to_fastq(bam, out_dir)

        assert result == out_dir / "extracted_reads.fq"
        assert result.read_text() == fastq_content

    @patch("muc_one_span.mapping.run_tool")
    def test_creates_output_dir(self, mock_run_tool, tmp_path):
        """bam_to_fastq creates the output directory if it does not exist."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "input.bam"
        bam.touch()
        out_dir = tmp_path / "nested" / "out"

        assert not out_dir.exists()
        bam_to_fastq(bam, out_dir)
        assert out_dir.exists()


class TestMapReads:
    """Tests for map_reads."""

    @patch("muc_one_span.mapping._run_mapping_pipeline")
    @patch("muc_one_span.mapping.run_tool_pipeline")
    def test_fastq_input_pipeline(self, mock_run_tool, mock_pipeline, tmp_path):
        """For FASTQ input, mapping pipeline and samtools index are called."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        out_dir = tmp_path / "out"

        map_reads(fastq, ref, out_dir, threads=2)

        # Pipeline should be called once (minimap2 | samtools sort)
        assert mock_pipeline.call_count == 1
        # samtools index should be called via run_tool
        index_cmd = mock_run_tool.call_args[0][0][0]
        assert index_cmd[:2] == ["samtools", "index"]

    @patch("muc_one_span.mapping._run_mapping_pipeline")
    @patch("muc_one_span.mapping.run_tool_pipeline")
    def test_bam_input_converts_to_fastq_first(self, mock_run_tool, mock_pipeline, tmp_path):
        """For BAM input, samtools fastq is called before the mapping pipeline."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "reads.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        out_dir = tmp_path / "out"

        map_reads(bam, ref, out_dir, threads=2)

        # First run_tool call should be samtools fastq (bam_to_fastq)
        first_cmd = mock_run_tool.call_args_list[0][0][0][0]
        assert first_cmd[:2] == ["samtools", "fastq"]

    @patch("muc_one_span.mapping._run_mapping_pipeline")
    @patch("muc_one_span.mapping.run_tool_pipeline")
    def test_returns_bam_path(self, mock_run_tool, mock_pipeline, tmp_path):
        """map_reads returns the sorted BAM path."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        result = map_reads(fastq, ref, tmp_path, threads=1)

        assert result == tmp_path / "mapping.bam"

    @patch("muc_one_span.mapping._run_mapping_pipeline")
    @patch("muc_one_span.mapping.run_tool_pipeline")
    def test_threads_passed_to_pipeline(self, mock_run_tool, mock_pipeline, tmp_path):
        """The threads parameter is passed to the mapping pipeline."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        map_reads(fastq, ref, tmp_path, threads=8)

        # Check threads were passed to _run_mapping_pipeline
        pipeline_call = mock_pipeline.call_args
        assert pipeline_call[0][3] == 8  # threads is the 4th positional arg

    @patch("muc_one_span.mapping._run_mapping_pipeline")
    @patch("muc_one_span.mapping.run_tool_pipeline")
    def test_no_intermediate_sam_file(self, mock_run_tool, mock_pipeline, tmp_path):
        """No intermediate SAM file is created (pipeline streams directly)."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        map_reads(fastq, ref, tmp_path, threads=1)
        assert not (tmp_path / "mapping.sam").exists()

    @patch("muc_one_span.mapping._run_mapping_pipeline")
    @patch("muc_one_span.mapping.run_tool_pipeline")
    def test_preset_passed_to_pipeline(self, mock_run_tool, mock_pipeline, tmp_path):
        """map_reads passes the preset parameter to _run_mapping_pipeline."""
        mock_run_tool.return_value = ""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        map_reads(fastq, ref, tmp_path, threads=2, preset="lr:hq")

        pipeline_call = mock_pipeline.call_args
        assert pipeline_call.kwargs["preset"] == "lr:hq"


class TestGetIdxstats:
    """Tests for get_idxstats."""

    @patch("muc_one_span.mapping.run_tool")
    def test_calls_samtools_idxstats(self, mock_run_tool, tmp_path):
        """get_idxstats calls samtools idxstats with the BAM path."""
        idxstats_output = "contig_60\t4120\t245\t0\n*\t0\t0\t50\n"
        mock_run_tool.return_value = idxstats_output
        bam = tmp_path / "mapping.bam"

        result = get_idxstats(bam)

        mock_run_tool.assert_called_once_with(["samtools", "idxstats", str(bam)])
        assert result == idxstats_output

    @patch("muc_one_span.mapping.run_tool")
    def test_returns_raw_output(self, mock_run_tool, tmp_path):
        """get_idxstats returns the raw string from run_tool unchanged."""
        raw = "contig_51\t3180\t115\t0\n"
        mock_run_tool.return_value = raw
        bam = tmp_path / "mapping.bam"

        assert get_idxstats(bam) == raw


@pytest.mark.parametrize("preset", ["map-hifi", "lr:hq", "map-ont"])
def test_pipeline_preserves_arguments_and_explicit_presets(mocker, preset):
    from muc_one_span.mapping import _run_mapping_pipeline

    run = mocker.patch("muc_one_span.mapping.run_tool_pipeline")
    _run_mapping_pipeline(
        Path("reads.fq"), Path("ref.fa"), Path("out.bam"), 3, preset=preset, timeout=27
    )
    assert run.call_args.args[0] == [
        ["minimap2", "-a", "-x", preset, "-t", "3", "ref.fa", "reads.fq"],
        ["samtools", "sort", "-@", "3", "-o", "out.bam"],
    ]
    assert run.call_args.kwargs["timeout"] == 27


def test_pipeline_default_preset_is_map_hifi(mocker):
    from muc_one_span.mapping import _run_mapping_pipeline

    run = mocker.patch("muc_one_span.mapping.run_tool_pipeline")
    _run_mapping_pipeline(Path("reads.fq"), Path("ref.fa"), Path("out.bam"), 1)
    assert run.call_args.args[0][0][3] == "map-hifi"
