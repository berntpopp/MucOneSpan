"""CLI pipeline execution and report fallback tests."""

from __future__ import annotations

from unittest.mock import patch

from click.testing import CliRunner

from open_pacmuci.cli import main


class TestRunSubcommand:
    """Tests for the run subcommand with all pipeline stages mocked."""

    def test_run_full_pipeline(self, tmp_path):
        """run subcommand completes successfully when all stages are mocked."""
        # Create a real input file (click checks exists=True)
        input_file = tmp_path / "reads.fastq"
        input_file.touch()
        output_dir = tmp_path / "results"

        # Create fake consensus FASTA files that will be read by the run function
        allele1_fa = tmp_path / "allele1.fa"
        allele1_fa.write_text(">allele_1_vntr\nACGTACGT\n")
        allele2_fa = tmp_path / "allele2.fa"
        allele2_fa.write_text(">allele_2_vntr\nACGTACGT\n")

        fake_alleles_result = {
            "homozygous": False,
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_51"],
            },
            "allele_2": {
                "length": 80,
                "reads": 150,
                "canonical_repeats": 71,
                "contig_name": "contig_71",
                "cluster_contigs": ["contig_71"],
            },
        }

        fake_classify_result = {
            "structure": "1-2-3-A-B-C-6-7",
            "repeats": ["A", "B", "C"],
            "mutations_detected": [],
            "allele_confidence": 0.95,
        }

        with (
            patch("open_pacmuci.tools.check_tools", return_value=True),
            patch(
                "open_pacmuci.tools.get_tool_versions",
                return_value={"minimap2": "2.24", "samtools": "1.17"},
            ),
            patch(
                "open_pacmuci.mapping.map_reads",
                return_value=tmp_path / "mapping.bam",
            ),
            patch(
                "open_pacmuci.mapping.get_idxstats",
                return_value="contig_51\t3120\t200\t0\ncontig_71\t4320\t150\t0\n*\t0\t0\t50\n",
            ),
            patch(
                "open_pacmuci.alleles.parse_idxstats",
                return_value={51: 200, 71: 150},
            ),
            patch(
                "open_pacmuci.alleles.detect_alleles",
                return_value=fake_alleles_result,
            ),
            patch(
                "open_pacmuci.calling.call_variants_per_allele",
                return_value={
                    "allele_1": tmp_path / "allele_1" / "variants.vcf.gz",
                    "allele_2": tmp_path / "allele_2" / "variants.vcf.gz",
                },
            ),
            patch(
                "open_pacmuci.consensus.build_consensus_per_allele",
                return_value={"allele_1": allele1_fa, "allele_2": allele2_fa},
            ),
            patch(
                "open_pacmuci.classify.classify_sequence",
                return_value=fake_classify_result,
            ),
            patch(
                "open_pacmuci.classify.validate_mutations_against_vcf",
                return_value=fake_classify_result,
            ),
            patch(
                "open_pacmuci.vcf.parse_vcf_variants",
                return_value=[],
            ),
        ):
            runner = CliRunner()
            result = runner.invoke(
                main,
                ["run", "--input", str(input_file), "--output-dir", str(output_dir)],
            )

        assert result.exit_code == 0, result.output
        assert "Pipeline complete" in result.output


class TestRunReportFallback:
    """Test --report flag graceful fallback when Jinja2 is missing."""

    def test_report_flag_warns_when_jinja2_missing(self, tmp_path):
        """--report warns and continues when Jinja2 is not installed."""
        input_file = tmp_path / "reads.fastq"
        input_file.touch()
        output_dir = tmp_path / "results"

        allele1_fa = tmp_path / "allele1.fa"
        allele1_fa.write_text(">a1\nACGT\n")

        fake_alleles = {
            "homozygous": False,
            "allele_1": {
                "length": 60,
                "reads": 100,
                "canonical_repeats": 51,
                "contig_name": "c51",
                "cluster_contigs": ["c51"],
            },
            "allele_2": {
                "length": 60,
                "reads": 100,
                "canonical_repeats": 51,
                "contig_name": "c51",
                "cluster_contigs": ["c51"],
            },
        }
        fake_cls = {
            "structure": "1 2 3 X 6 7 8 9",
            "repeats": [],
            "mutations_detected": [],
            "allele_confidence": 1.0,
        }

        with (
            patch("open_pacmuci.tools.check_tools"),
            patch("open_pacmuci.tools.get_tool_versions", return_value={}),
            patch("open_pacmuci.mapping.map_reads", return_value=tmp_path / "m.bam"),
            patch("open_pacmuci.mapping.get_idxstats", return_value="c51\t3060\t100\t0\n"),
            patch("open_pacmuci.alleles.parse_idxstats", return_value={51: 100}),
            patch("open_pacmuci.alleles.detect_alleles", return_value=fake_alleles),
            patch(
                "open_pacmuci.calling.call_variants_per_allele",
                return_value={"allele_1": tmp_path / "a.vcf.gz"},
            ),
            patch(
                "open_pacmuci.consensus.build_consensus_per_allele",
                return_value={"allele_1": allele1_fa},
            ),
            patch("open_pacmuci.classify.classify_sequence", return_value=fake_cls),
            patch("open_pacmuci.classify.validate_mutations_against_vcf", return_value=fake_cls),
            patch("open_pacmuci.vcf.parse_vcf_variants", return_value=[]),
            patch("open_pacmuci.cli._bundled_reference", return_value=tmp_path / "ref.fa"),
            patch.dict("sys.modules", {"jinja2": None}),
            patch("open_pacmuci.report._HAS_JINJA2", False),
        ):
            runner = CliRunner()
            result = runner.invoke(
                main,
                ["run", "--input", str(input_file), "--output-dir", str(output_dir), "--report"],
            )

        assert result.exit_code == 0, result.output
        assert "Pipeline complete" in result.output
