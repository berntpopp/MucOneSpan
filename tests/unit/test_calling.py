"""Mocked unit tests for calling module (Clair3 / bcftools wrappers)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from muc_one_span.calling import (
    _extract_and_remap_reads,
    call_variants_per_allele,
    disambiguate_same_length_alleles,
    extract_allele_reads,
    run_clair3,
)


class TestExtractAlleleReads:
    """Tests for extract_allele_reads."""

    @patch("muc_one_span.calling.run_tool")
    def test_single_contig_string(self, mock_run_tool, tmp_path):
        """Accepts a single contig name as a string and calls samtools view."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()

        extract_allele_reads(bam, "contig_51", tmp_path / "out")

        view_call = mock_run_tool.call_args_list[0][0][0]
        assert view_call[:2] == ["samtools", "view"]
        assert "contig_51" in view_call

    @patch("muc_one_span.calling.run_tool")
    def test_list_of_contigs(self, mock_run_tool, tmp_path):
        """Accepts a list of contig names and passes them all to samtools view."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        contigs = ["contig_50", "contig_51", "contig_52"]

        extract_allele_reads(bam, contigs, tmp_path / "out")

        view_call = mock_run_tool.call_args_list[0][0][0]
        for c in contigs:
            assert c in view_call

    @patch("muc_one_span.calling.run_tool")
    def test_returns_allele_bam_path(self, mock_run_tool, tmp_path):
        """Returns path to allele_reads.bam inside output dir."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        out_dir = tmp_path / "out"

        result = extract_allele_reads(bam, "contig_51", out_dir)

        assert result == out_dir / "allele_reads.bam"

    @patch("muc_one_span.calling.run_tool")
    def test_indexes_result_bam(self, mock_run_tool, tmp_path):
        """samtools index is called on the output BAM."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        out_dir = tmp_path / "out"

        extract_allele_reads(bam, "contig_51", out_dir)

        index_call = mock_run_tool.call_args_list[1][0][0]
        assert index_call[:2] == ["samtools", "index"]

    @patch("muc_one_span.calling.run_tool")
    def test_creates_output_dir(self, mock_run_tool, tmp_path):
        """Creates the output directory if it does not exist."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        out_dir = tmp_path / "nested" / "out"

        assert not out_dir.exists()
        extract_allele_reads(bam, "contig_51", out_dir)
        assert out_dir.exists()


class TestRunClair3:
    """Tests for run_clair3."""

    @patch("muc_one_span.calling.run_tool")
    def test_builds_correct_command(self, mock_run_tool, tmp_path):
        """run_clair3 passes required flags to run_clair3.sh."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        run_clair3(bam, ref, out_dir, platform="hifi", threads=4)

        cmd = mock_run_tool.call_args[0][0]
        assert cmd[0] == "run_clair3.sh"
        assert f"--bam_fn={bam}" in cmd
        assert f"--ref_fn={ref}" in cmd
        assert "--include_all_ctgs" in cmd
        assert "--platform=hifi" in cmd
        assert "--threads=4" in cmd

    @patch("muc_one_span.calling.run_tool")
    def test_model_path_appended_when_given(self, mock_run_tool, tmp_path):
        """--model_path flag is added when model_path is non-empty."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        run_clair3(bam, ref, out_dir, model_path="/models/hifi")

        cmd = mock_run_tool.call_args[0][0]
        assert "--model_path=/models/hifi" in cmd

    @patch("muc_one_span.calling.run_tool")
    def test_model_path_omitted_when_empty(self, mock_run_tool, tmp_path):
        """--model_path is not added when model_path is empty string."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        run_clair3(bam, ref, out_dir, model_path="")

        cmd = mock_run_tool.call_args[0][0]
        assert not any(a.startswith("--model_path") for a in cmd)

    @patch("muc_one_span.calling.run_tool")
    def test_returns_vcf_path(self, mock_run_tool, tmp_path):
        """run_clair3 returns path to merge_output.vcf.gz."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3"

        result = run_clair3(bam, ref, out_dir)

        assert result == out_dir / "merge_output.vcf.gz"

    @patch("muc_one_span.calling.run_tool")
    def test_creates_output_dir(self, mock_run_tool, tmp_path):
        """run_clair3 creates the output directory if it does not exist."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "allele.bam"
        ref = tmp_path / "ref.fa"
        out_dir = tmp_path / "clair3" / "nested"

        assert not out_dir.exists()
        run_clair3(bam, ref, out_dir)
        assert out_dir.exists()


class TestCallVariantsPerAllele:
    """Tests for call_variants_per_allele."""

    def _make_alleles_heterozygous(self):
        return {
            "homozygous": False,
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_50", "contig_51", "contig_52"],
            },
            "allele_2": {
                "length": 80,
                "reads": 150,
                "canonical_repeats": 71,
                "contig_name": "contig_71",
                "cluster_contigs": ["contig_70", "contig_71", "contig_72"],
            },
        }

    def _make_alleles_homozygous(self):
        return {
            "homozygous": True,
            "allele_1": {
                "length": 60,
                "reads": 400,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_50", "contig_51", "contig_52"],
            },
            "allele_2": {
                "length": 60,
                "reads": 0,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_50", "contig_51", "contig_52"],
            },
        }

    @patch("muc_one_span.vcf.run_tool")
    @patch("muc_one_span.calling.run_tool")
    def test_heterozygous_processes_both_alleles(self, mock_run_tool, mock_vcf_tool, tmp_path):
        """For a heterozygous sample, both allele_1 and allele_2 are processed."""
        mock_run_tool.return_value = ""
        mock_vcf_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        alleles = self._make_alleles_heterozygous()

        result = call_variants_per_allele(bam, ref, alleles, tmp_path)

        assert "allele_1" in result
        assert "allele_2" in result

    @patch("muc_one_span.vcf.run_tool")
    @patch("muc_one_span.calling.run_tool")
    def test_homozygous_skips_allele_2(self, mock_run_tool, mock_vcf_tool, tmp_path):
        """For a homozygous sample, allele_2 is skipped."""
        mock_run_tool.return_value = ""
        mock_vcf_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        alleles = self._make_alleles_homozygous()

        result = call_variants_per_allele(bam, ref, alleles, tmp_path)

        assert "allele_1" in result
        assert "allele_2" not in result

    @patch("muc_one_span.vcf.run_tool")
    @patch("muc_one_span.calling.run_tool")
    def test_returns_dict_of_vcf_paths(self, mock_run_tool, mock_vcf_tool, tmp_path):
        """Results map allele keys to Path objects."""
        mock_run_tool.return_value = ""
        mock_vcf_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        alleles = self._make_alleles_heterozygous()

        result = call_variants_per_allele(bam, ref, alleles, tmp_path)

        for _key, path in result.items():
            assert isinstance(path, Path)

    @patch("muc_one_span.vcf.run_tool")
    @patch("muc_one_span.calling.run_tool")
    def test_fallback_contig_name_from_length(self, mock_run_tool, mock_vcf_tool, tmp_path):
        """When contig_name is absent, falls back to contig_<length>."""
        mock_run_tool.return_value = ""
        mock_vcf_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        # Alleles dict without explicit contig_name (legacy format)
        alleles = {
            "homozygous": True,
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                "cluster_contigs": ["contig_51"],
                # no "contig_name" key
            },
            "allele_2": {
                "length": 60,
                "reads": 0,
                "canonical_repeats": 51,
                "cluster_contigs": ["contig_51"],
            },
        }

        # Should not raise
        result = call_variants_per_allele(bam, ref, alleles, tmp_path)
        assert "allele_1" in result

    @patch("muc_one_span.vcf.run_tool")
    @patch("muc_one_span.calling.run_tool")
    def test_platform_and_preset_threaded_through(self, mock_run_tool, mock_vcf_tool, tmp_path):
        """platform and preset are forwarded: minimap2 gets -x lr:hq, clair3 gets --platform=ont."""
        mock_run_tool.return_value = ""
        mock_vcf_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        alleles = self._make_alleles_heterozygous()

        call_variants_per_allele(bam, ref, alleles, tmp_path, platform="ont", preset="lr:hq")

        all_calls = [c[0][0] for c in mock_run_tool.call_args_list]
        minimap2_calls = [cmd for cmd in all_calls if cmd[0] == "minimap2"]
        clair3_calls = [cmd for cmd in all_calls if cmd[0] == "run_clair3.sh"]

        assert minimap2_calls, "minimap2 was not called"
        minimap2_cmd = minimap2_calls[0]
        x_idx = minimap2_cmd.index("-x")
        assert minimap2_cmd[x_idx + 1] == "lr:hq"

        assert clair3_calls, "run_clair3.sh was not called"
        assert any("--platform=ont" in cmd for cmd in clair3_calls)


class TestExtractAndRemapReads:
    """Tests for the private _extract_and_remap_reads helper."""

    @patch("muc_one_span.calling.run_tool")
    def test_remaps_to_peak_contig(self, mock_run_tool, tmp_path):
        """Remapping pipeline is: extract→fastq→faidx (extract contig)→faidx (index)→minimap2→sort→index."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        _extract_and_remap_reads(
            bam,
            ["contig_50", "contig_51", "contig_52"],
            "contig_51",
            ref,
            tmp_path / "out",
            threads=2,
        )

        tool_names = [c[0][0][0] for c in mock_run_tool.call_args_list]
        assert "minimap2" in tool_names
        # samtools sort should be present
        sort_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["samtools", "sort"]
        ]
        assert sort_calls

    @patch("muc_one_span.calling.run_tool")
    def test_creates_contig_fasta_for_reference(self, mock_run_tool, tmp_path):
        """faidx is called to extract the peak contig as a mini-reference."""
        mock_run_tool.return_value = ">contig_51\nACGT\n"
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        out_dir = tmp_path / "out"

        _extract_and_remap_reads(bam, ["contig_51"], "contig_51", ref, out_dir, threads=1)

        faidx_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["samtools", "faidx"]
        ]
        # At least one faidx call with the peak contig name
        contig_faidx = [c for c in faidx_calls if "contig_51" in c]
        assert contig_faidx

    @patch("muc_one_span.calling.run_tool")
    def test_preset_passed_to_minimap2(self, mock_run_tool, tmp_path):
        """When preset='lr:hq' is given, minimap2 is called with -x lr:hq."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        _extract_and_remap_reads(
            bam,
            ["contig_51"],
            "contig_51",
            ref,
            tmp_path / "out",
            threads=2,
            preset="lr:hq",
        )

        all_calls = [c[0][0] for c in mock_run_tool.call_args_list]
        minimap2_calls = [cmd for cmd in all_calls if cmd[0] == "minimap2"]
        assert minimap2_calls, "minimap2 was not called"
        minimap2_cmd = minimap2_calls[0]
        x_idx = minimap2_cmd.index("-x")
        assert minimap2_cmd[x_idx + 1] == "lr:hq"

    @patch("muc_one_span.calling.run_tool")
    def test_preset_defaults_to_map_hifi(self, mock_run_tool, tmp_path):
        """When preset is not given, minimap2 is called with -x map-hifi."""
        mock_run_tool.return_value = ""
        bam = tmp_path / "mapping.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        _extract_and_remap_reads(
            bam,
            ["contig_51"],
            "contig_51",
            ref,
            tmp_path / "out",
            threads=2,
        )

        all_calls = [c[0][0] for c in mock_run_tool.call_args_list]
        minimap2_calls = [cmd for cmd in all_calls if cmd[0] == "minimap2"]
        assert minimap2_calls, "minimap2 was not called"
        minimap2_cmd = minimap2_calls[0]
        x_idx = minimap2_cmd.index("-x")
        assert minimap2_cmd[x_idx + 1] == "map-hifi"


class TestDisambiguateSameLengthAlleles:
    """Tests for disambiguate_same_length_alleles."""

    @patch("muc_one_span.vcf.run_tool", return_value="")
    @patch("muc_one_span.calling.run_tool", return_value="")
    @patch("muc_one_span.calling.parse_vcf_genotypes", return_value=[])
    def test_no_het_variants_returns_homozygous(self, mock_geno, mock_run, mock_vcf_tool, tmp_path):
        alleles = {
            "allele_1": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
            "allele_2": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
        }
        result = disambiguate_same_length_alleles(
            tmp_path / "bam", tmp_path / "ref.fa", alleles, tmp_path
        )
        assert "allele_1" in result
        assert result.get("homozygous") is True

    @patch("muc_one_span.vcf.run_tool", return_value="")
    @patch("muc_one_span.calling.run_tool", return_value="")
    @patch(
        "muc_one_span.calling.parse_vcf_genotypes",
        return_value=[{"chrom": "c", "pos": 100, "ref": "A", "alt": "T", "genotype": "0/1"}],
    )
    def test_het_variants_returns_two_alleles(self, mock_geno, mock_run, mock_vcf_tool, tmp_path):
        alleles = {
            "allele_1": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
            "allele_2": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
        }
        result = disambiguate_same_length_alleles(
            tmp_path / "bam", tmp_path / "ref.fa", alleles, tmp_path
        )
        assert "allele_1" in result
        assert "allele_2" in result
        assert result.get("homozygous") is False

    @patch("muc_one_span.vcf.run_tool", return_value="")
    @patch("muc_one_span.calling.run_tool", return_value="")
    @patch("muc_one_span.calling.parse_vcf_genotypes", return_value=[])
    def test_platform_and_preset_threaded_through(
        self, mock_geno, mock_run, mock_vcf_tool, tmp_path
    ):
        """platform and preset are forwarded: minimap2 gets -x lr:hq, clair3 gets --platform=ont."""
        alleles = {
            "allele_1": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
            "allele_2": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
        }
        disambiguate_same_length_alleles(
            tmp_path / "bam",
            tmp_path / "ref.fa",
            alleles,
            tmp_path,
            platform="ont",
            preset="lr:hq",
        )

        all_calls = [c[0][0] for c in mock_run.call_args_list]
        minimap2_calls = [cmd for cmd in all_calls if cmd[0] == "minimap2"]
        clair3_calls = [cmd for cmd in all_calls if cmd[0] == "run_clair3.sh"]

        assert minimap2_calls, "minimap2 was not called"
        minimap2_cmd = minimap2_calls[0]
        x_idx = minimap2_cmd.index("-x")
        assert minimap2_cmd[x_idx + 1] == "lr:hq"

        assert clair3_calls, "run_clair3.sh was not called"
        assert any("--platform=ont" in cmd for cmd in clair3_calls)
