"""Unit tests for vcf module (filter_vcf, parse_vcf_genotypes, parse_vcf_variants)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.vcf import filter_vcf, parse_vcf_genotypes, parse_vcf_variants


class TestFilterVcf:
    """Tests for filter_vcf."""

    @patch("muc_one_span.vcf.run_tool")
    def test_calls_bcftools_norm_then_view_then_index(self, mock_run_tool, tmp_path):
        """filter_vcf runs bcftools norm, bcftools view, then bcftools index."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        filter_vcf(vcf, ref, out_dir)

        calls = mock_run_tool.call_args_list
        assert calls[0][0][0][:2] == ["bcftools", "norm"]
        assert calls[1][0][0][:2] == ["bcftools", "view"]
        assert calls[2][0][0][:2] == ["bcftools", "index"]

    @patch("muc_one_span.vcf.run_tool")
    def test_norm_uses_reference(self, mock_run_tool, tmp_path):
        """bcftools norm -f flag receives the reference path."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        filter_vcf(vcf, ref, out_dir)

        norm_cmd = mock_run_tool.call_args_list[0][0][0]
        assert "-f" in norm_cmd
        ref_idx = norm_cmd.index("-f")
        assert norm_cmd[ref_idx + 1] == str(ref)

    @patch("muc_one_span.vcf.run_tool")
    def test_view_filters_pass(self, mock_run_tool, tmp_path):
        """bcftools view uses -f PASS to keep only passing variants."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        filter_vcf(vcf, ref, out_dir)

        view_cmd = mock_run_tool.call_args_list[1][0][0]
        assert "-f" in view_cmd
        f_idx = view_cmd.index("-f")
        assert view_cmd[f_idx + 1] == "PASS"

    @patch("muc_one_span.vcf.run_tool")
    def test_returns_variants_vcf_path(self, mock_run_tool, tmp_path):
        """filter_vcf returns path to variants.vcf.gz."""
        mock_run_tool.return_value = ""
        vcf = tmp_path / "raw.vcf.gz"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "filtered"

        result = filter_vcf(vcf, ref, out_dir)

        assert result == out_dir / "variants.vcf.gz"


class TestFilterVcfEmptyVcf:
    """Tests for filter_vcf handling empty VCFs (issue #8)."""

    @patch("muc_one_span.vcf.run_tool")
    def test_empty_vcf_skips_quality_filter(self, mock_run_tool, tmp_path):
        """Empty VCF (no records) skips -i filter to avoid INFO/DP crash."""
        # bcftools norm returns empty output, then bcftools view should NOT
        # include -i filter since there are no records to filter
        norm_vcf = tmp_path / "normalized.vcf.gz"

        def side_effect(cmd):
            # After norm runs, create the normalized file with header only
            if cmd[:2] == ["bcftools", "norm"]:
                norm_vcf.write_bytes(b"")  # empty file
            return ""

        mock_run_tool.side_effect = side_effect
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        filter_vcf(vcf, ref, tmp_path, min_qual=15.0, min_dp=5)

        # bcftools view should NOT have -i flag (empty VCF)
        view_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["bcftools", "view"]
        ]
        assert len(view_calls) == 1
        assert "-i" not in view_calls[0]

    @patch("muc_one_span.vcf.run_tool")
    def test_nonempty_vcf_uses_qual_only_filter(self, mock_run_tool, tmp_path):
        """Non-empty VCF uses QUAL filter only (no INFO/DP which may not exist)."""
        norm_vcf = tmp_path / "normalized.vcf.gz"

        def side_effect(cmd):
            if cmd[:2] == ["bcftools", "norm"]:
                # Write a non-empty file to simulate records
                norm_vcf.write_bytes(b"\x00" * 100)
            return ""

        mock_run_tool.side_effect = side_effect
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        filter_vcf(vcf, ref, tmp_path, min_qual=15.0, min_dp=5)

        view_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["bcftools", "view"]
        ]
        assert len(view_calls) == 1
        view_cmd = view_calls[0]
        assert "-i" in view_cmd
        i_idx = view_cmd.index("-i")
        expr = view_cmd[i_idx + 1]
        assert "QUAL" in expr
        # Should NOT reference INFO/DP (may not exist in Clair3 output)
        assert "INFO/DP" not in expr


class TestFilterVcfQuality:
    """Tests for VCF quality filter parameters."""

    @patch("muc_one_span.vcf.run_tool")
    def test_filter_vcf_includes_quality_expression(self, mock_run_tool, tmp_path):
        """filter_vcf passes QUAL filter to bcftools view for non-empty VCFs."""
        norm_vcf = tmp_path / "normalized.vcf.gz"

        def side_effect(cmd):
            if cmd[:2] == ["bcftools", "norm"]:
                norm_vcf.write_bytes(b"\x00" * 100)  # non-empty
            return ""

        mock_run_tool.side_effect = side_effect
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        filter_vcf(vcf, ref, tmp_path, min_qual=15.0, min_dp=5)

        # Find the bcftools view call
        view_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["bcftools", "view"]
        ]
        assert len(view_calls) == 1
        view_cmd = view_calls[0]
        assert "-i" in view_cmd
        i_idx = view_cmd.index("-i")
        expr = view_cmd[i_idx + 1]
        assert "QUAL" in expr

    @patch("muc_one_span.vcf.run_tool", return_value="")
    def test_filter_vcf_default_params(self, mock_run_tool, tmp_path):
        """filter_vcf works with default parameters (backward compatible)."""
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        # Should not raise with no extra args
        filter_vcf(vcf, ref, tmp_path)


@pytest.mark.parametrize("parser", [parse_vcf_genotypes, parse_vcf_variants])
def test_query_error_is_not_an_empty_call(parser):
    with (
        patch("muc_one_span.vcf.run_tool", side_effect=RuntimeError("query failed")),
        pytest.raises(RuntimeError, match="query failed"),
    ):
        parser(Path("fake.vcf"))


@pytest.mark.parametrize("parser", [parse_vcf_genotypes, parse_vcf_variants])
def test_retains_variant_identity_missing_quality_and_phase(parser):
    with patch(
        "muc_one_span.vcf.run_tool", side_effect=["SAMPLE\n", "chr1\t2\tA\tC,G\t.\t1|2\t17\n"]
    ):
        assert parser(Path("fake.vcf")) == [
            {
                "chrom": "chr1",
                "pos": 2,
                "ref": "A",
                "alt": "C,G",
                "qual": None,
                "genotype": "1|2",
                "phase_set": "17",
                "sample": "SAMPLE",
            }
        ]


@pytest.mark.parametrize(
    "record",
    [
        "bad",
        "c\t0\tA\tC\t3\t0/1\t.",
        "c\t2\tA\tC\tnan\t0/1\t.",
        "c\t2\tA\tC\t3\t0/x\t.",
        "c\t2\tA\tC\t3\t0/2\t.",
    ],
)
def test_malformed_query_fails_visibly(record):
    with (
        patch("muc_one_span.vcf.run_tool", side_effect=["SAMPLE\n", record]),
        pytest.raises(ValueError),
    ):
        parse_vcf_variants(Path("fake.vcf"))


def test_successful_empty_vcf_retains_empty_result():
    with patch("muc_one_span.vcf.run_tool", side_effect=["SAMPLE\n", ""]):
        assert parse_vcf_variants(Path("fake.vcf")) == []


def test_multiple_samples_require_selection():
    with (
        patch("muc_one_span.vcf.run_tool", return_value="A\nB\n"),
        pytest.raises(ValueError, match="sample"),
    ):
        parse_vcf_variants(Path("fake.vcf"))


def test_explicit_sample_is_selected():
    with patch("muc_one_span.vcf.run_tool", side_effect=["A\nB\n", ""]) as run:
        assert parse_vcf_variants(Path("fake.vcf"), sample="B") == []
        assert run.call_args.args[0][-3:] == ["-s", "B", "fake.vcf"]


def test_filter_vcf_haploid_majority(tmp_path):
    vcf = tmp_path / "input.vcf.gz"
    vcf.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()
    out_dir = tmp_path / "out"

    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "contig_1\t10\t.\tC\tG\t15.0\tPASS\t.\tGT:AF\t0/1:0.8\n"
        "contig_1\t20\t.\tA\tT\t12.0\tPASS\t.\tGT:AF\t0/1:0.1\n"
        "contig_1\t30\t.\tG\tC\t14.0\tPASS\t.\tGT:AF\t0/1:0.35\n"
    )

    written_lines: list[str] = []

    def side_effect(cmd):
        if cmd[0] == "bcftools" and cmd[1] == "norm":
            (out_dir / "normalized.vcf.gz").write_bytes(b"data")
            return ""
        if cmd[0] == "bcftools" and cmd[1] == "view" and str(cmd[-1]).endswith("normalized.vcf.gz"):
            (out_dir / "variants.vcf.gz").write_bytes(b"data")
            return ""
        if cmd[0] == "bcftools" and cmd[1] == "view" and str(cmd[-1]).endswith("variants.vcf.gz"):
            return vcf_content
        if cmd[0] == "bcftools" and cmd[1] == "view" and str(cmd[-1]).endswith("mod_haploid.vcf"):
            mod_file = Path(cmd[-1])
            written_lines.extend(mod_file.read_text().splitlines())
            (out_dir / "variants.vcf.gz").touch()
            return ""
        return ""

    with patch("muc_one_span.vcf.run_tool", side_effect=side_effect):
        filter_vcf(vcf, ref, out_dir, min_qual=5.0, haploid_majority=True)

    assert any("0/1:0.8" not in line and "1/1:0.8" in line for line in written_lines)
    assert any("0/1:0.1" not in line and "0/0:0.1" in line for line in written_lines)
    assert any("0/1:0.35" in line for line in written_lines)
