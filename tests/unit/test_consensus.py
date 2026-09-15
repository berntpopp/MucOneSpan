"""Mocked unit tests for consensus module (bcftools wrappers)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.consensus import build_consensus, build_consensus_per_allele, trim_flanking


@pytest.fixture(autouse=True)
def selected_sample():
    # Sample header inspection is an external bcftools operation.
    with patch("muc_one_span.consensus.select_vcf_sample", return_value="SAMPLE"):
        yield


class TestBuildConsensus:
    """Tests for build_consensus."""

    @patch("muc_one_span.consensus.run_tool")
    def test_calls_bcftools_consensus(self, mock_run_tool, tmp_path):
        """build_consensus calls bcftools consensus with correct arguments."""
        fasta_content = ">contig_51\nACGTACGT\n"
        mock_run_tool.return_value = fasta_content
        ref = tmp_path / "ref.fa"
        ref.touch()
        vcf = tmp_path / "variants.vcf.gz"
        vcf.touch()
        output = tmp_path / "consensus.fa"

        build_consensus(ref, vcf, output)

        cmd = mock_run_tool.call_args[0][0]
        assert cmd[:2] == ["bcftools", "consensus"]
        assert str(ref) in cmd
        assert str(vcf) in cmd

    @patch("muc_one_span.consensus.run_tool")
    def test_writes_stdout_to_output_path(self, mock_run_tool, tmp_path):
        """build_consensus writes the bcftools stdout to the output file."""
        fasta_content = ">contig_51\nACGTACGT\n"
        mock_run_tool.return_value = fasta_content
        ref = tmp_path / "ref.fa"
        ref.touch()
        vcf = tmp_path / "variants.vcf.gz"
        vcf.touch()
        output = tmp_path / "consensus.fa"

        result = build_consensus(ref, vcf, output)

        assert output.read_text() == fasta_content
        assert result == output

    @patch("muc_one_span.consensus.run_tool")
    def test_creates_parent_directory(self, mock_run_tool, tmp_path):
        """build_consensus creates the parent directory of output_path if needed."""
        mock_run_tool.return_value = ">x\nACGT\n"
        ref = tmp_path / "ref.fa"
        ref.touch()
        vcf = tmp_path / "variants.vcf.gz"
        vcf.touch()
        output = tmp_path / "nested" / "dir" / "consensus.fa"

        assert not output.parent.exists()
        build_consensus(ref, vcf, output)
        assert output.parent.exists()

    @patch("muc_one_span.consensus.run_tool")
    def test_reference_passed_with_f_flag(self, mock_run_tool, tmp_path):
        """bcftools consensus uses -f flag for the reference."""
        mock_run_tool.return_value = ">x\nACGT\n"
        ref = tmp_path / "ref.fa"
        ref.touch()
        vcf = tmp_path / "variants.vcf.gz"
        vcf.touch()
        output = tmp_path / "out.fa"

        build_consensus(ref, vcf, output)

        cmd = mock_run_tool.call_args[0][0]
        f_idx = cmd.index("-f")
        assert cmd[f_idx + 1] == str(ref)


class TestTrimFlanking:
    """Tests for trim_flanking (pure Python, no mocking needed)."""

    def test_removes_flanks_from_both_ends(self, tmp_path):
        """Correct flanking bp are trimmed from both ends of the sequence."""
        sequence = "A" * 100 + "VNTR_SEQ" + "A" * 100
        fasta = tmp_path / "full.fa"
        fasta.write_text(f">contig_51\n{sequence}\n")
        output = tmp_path / "trimmed.fa"

        trim_flanking(fasta, 100, output)

        lines = output.read_text().strip().splitlines()
        trimmed_seq = "".join(ln for ln in lines if not ln.startswith(">"))
        assert trimmed_seq == "VNTR_SEQ"

    def test_output_header_has_vntr_suffix(self, tmp_path):
        """Output FASTA header gets a _vntr suffix."""
        fasta = tmp_path / "full.fa"
        fasta.write_text(">contig_51\n" + "A" * 10 + "CORE" + "A" * 10 + "\n")
        output = tmp_path / "trimmed.fa"

        trim_flanking(fasta, 10, output)

        header = output.read_text().strip().splitlines()[0]
        assert header.endswith("_vntr")

    def test_zero_flank_length_keeps_full_sequence(self, tmp_path):
        """With flank_length=0, the entire sequence is kept."""
        sequence = "ACGTACGT"
        fasta = tmp_path / "full.fa"
        fasta.write_text(f">contig_51\n{sequence}\n")
        output = tmp_path / "trimmed.fa"

        trim_flanking(fasta, 0, output)

        lines = output.read_text().strip().splitlines()
        trimmed_seq = "".join(ln for ln in lines if not ln.startswith(">"))
        assert trimmed_seq == sequence

    def test_returns_output_path(self, tmp_path):
        """trim_flanking returns the output path."""
        fasta = tmp_path / "full.fa"
        fasta.write_text(">c\n" + "A" * 20 + "\n")
        output = tmp_path / "out.fa"

        result = trim_flanking(fasta, 5, output)

        assert result == output

    def test_handles_fasta_without_header(self, tmp_path):
        """Falls back to >consensus when no header line is present."""
        fasta = tmp_path / "full.fa"
        fasta.write_text("ACGT")
        output = tmp_path / "out.fa"

        trim_flanking(fasta, 0, output)

        header = output.read_text().strip().splitlines()[0]
        assert header == ">consensus_vntr"


class TestAnchorBasedTrim:
    """Tests for anchor-based flanking trim."""

    def test_anchor_trim_handles_indel_in_flank(self, tmp_path):
        """Anchor trim finds correct boundary despite indel in flanking."""
        from muc_one_span.config import load_repeat_dictionary

        rd = load_repeat_dictionary()
        # Build a sequence: left_flank (with 1bp insertion) + pre-repeat 1 + X + after-repeat 9 + right_flank
        left_flank = rd.flanking_left[-500:]
        right_flank = rd.flanking_right[:500]
        vntr = rd.repeats["1"] + rd.repeats["X"] + rd.repeats["9"]
        # Insert a fake base in the left flank (simulating Clair3 false positive)
        corrupted_left = left_flank[:250] + "A" + left_flank[250:]
        sequence = corrupted_left + vntr + right_flank  # 501 + vntr + 500

        fasta = tmp_path / "full.fa"
        fasta.write_text(f">contig\n{sequence}\n")
        output = tmp_path / "trimmed.fa"

        trim_flanking(fasta, 500, output, repeat_dict=rd)

        trimmed = output.read_text().strip().splitlines()
        trimmed_seq = "".join(ln for ln in trimmed if not ln.startswith(">"))
        # Should contain the start of repeat "1"
        assert rd.repeats["1"][:30] in trimmed_seq

    def test_anchor_trim_fallback_without_repeat_dict(self, tmp_path):
        """Without repeat_dict, falls back to fixed-position trim."""
        sequence = "A" * 100 + "VNTR_SEQ" + "A" * 100
        fasta = tmp_path / "full.fa"
        fasta.write_text(f">contig\n{sequence}\n")
        output = tmp_path / "trimmed.fa"

        trim_flanking(fasta, 100, output, repeat_dict=None)

        trimmed = output.read_text().strip().splitlines()
        trimmed_seq = "".join(ln for ln in trimmed if not ln.startswith(">"))
        assert trimmed_seq == "VNTR_SEQ"


class TestBuildConsensusPerAllele:
    """Tests for build_consensus_per_allele."""

    def _alleles(self):
        return {
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

    @patch("muc_one_span.consensus.run_tool")
    def test_processes_all_allele_keys_in_vcf_paths(self, mock_run_tool, tmp_path):
        """build_consensus_per_allele handles every key in vcf_paths."""
        flanks = 50
        # Fake consensus output: 2*flanks + some core
        core = "GCGCGCGC"
        fake_fasta = f">contig_51\n{'A' * flanks}{core}{'A' * flanks}\n"
        mock_run_tool.return_value = fake_fasta

        vcf1 = tmp_path / "variants1.vcf.gz"
        vcf1.touch()
        vcf2 = tmp_path / "variants2.vcf.gz"
        vcf2.touch()
        vcf_paths = {"allele_1": vcf1, "allele_2": vcf2}
        alleles = self._alleles()
        out_dir = tmp_path / "consensus"

        result = build_consensus_per_allele(
            tmp_path / "ref.fa", vcf_paths, alleles, out_dir, flank_length=flanks
        )

        assert "allele_1" in result
        assert "allele_2" in result

    @patch("muc_one_span.consensus.run_tool")
    def test_calls_samtools_faidx_for_each_allele(self, mock_run_tool, tmp_path):
        """samtools faidx is called to extract each allele's contig."""
        flanks = 10
        mock_run_tool.return_value = f">contig_51\n{'A' * (flanks * 2 + 8)}\n"

        vcf1 = tmp_path / "v.vcf.gz"
        vcf1.touch()
        vcf_paths = {"allele_1": vcf1}
        alleles = {
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_51"],
            }
        }

        build_consensus_per_allele(
            tmp_path / "ref.fa", vcf_paths, alleles, tmp_path / "out", flank_length=flanks
        )

        faidx_calls = [
            c[0][0] for c in mock_run_tool.call_args_list if c[0][0][:2] == ["samtools", "faidx"]
        ]
        # At least one faidx extracting the contig
        assert any("contig_51" in c for c in faidx_calls)

    @patch("muc_one_span.consensus.run_tool")
    def test_returns_paths_to_trimmed_fastas(self, mock_run_tool, tmp_path):
        """Returns paths to the trimmed (VNTR-only) FASTA files."""
        flanks = 10
        mock_run_tool.return_value = f">c\n{'A' * (flanks * 2 + 4)}\n"

        vcf1 = tmp_path / "v.vcf.gz"
        vcf1.touch()
        alleles = {
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_51"],
            }
        }

        result = build_consensus_per_allele(
            tmp_path / "ref.fa",
            {"allele_1": vcf1},
            alleles,
            tmp_path / "out",
            flank_length=flanks,
        )

        assert isinstance(result["allele_1"], Path)
        # Should be the trimmed file, not the full consensus
        assert "full" not in result["allele_1"].name

    @patch("muc_one_span.consensus.run_tool")
    def test_fallback_contig_name_from_length(self, mock_run_tool, tmp_path):
        """Falls back to contig_<length> when contig_name key is absent."""
        flanks = 5
        mock_run_tool.return_value = f">c\n{'A' * (flanks * 2 + 4)}\n"

        vcf1 = tmp_path / "v.vcf.gz"
        vcf1.touch()
        alleles = {
            "allele_1": {
                "length": 60,
                "reads": 200,
                "canonical_repeats": 51,
                # no "contig_name" key
                "cluster_contigs": ["contig_51"],
            }
        }

        # Should not raise
        result = build_consensus_per_allele(
            tmp_path / "ref.fa",
            {"allele_1": vcf1},
            alleles,
            tmp_path / "out",
            flank_length=flanks,
        )
        assert "allele_1" in result


def test_consensus_uses_explicit_sample_and_genotype_haplotype(tmp_path):
    with patch("muc_one_span.consensus.run_tool", return_value=">c\nAC\n") as run:
        build_consensus(
            Path("ref.fa"), Path("calls.vcf"), tmp_path / "out.fa", sample="SAMPLE", haplotype=2
        )
        command = run.call_args.args[0]
        assert command[command.index("-s") + 1] == "SAMPLE"
        assert command[command.index("-H") + 1] == "2"


def test_trim_records_actual_anchor_interval(tmp_path):
    from muc_one_span.config import load_repeat_dictionary

    rd = load_repeat_dictionary()
    left = rd.flanking_left[-500:]
    vntr = rd.repeats["1"] + rd.repeats["X"] + rd.repeats["9"]
    sequence = left[:250] + "A" + left[250:] + vntr + rd.flanking_right[:500]
    full = tmp_path / "full.fa"
    full.write_text(f">c\n{sequence}\n")
    context = {}
    trim_flanking(full, 500, tmp_path / "trim.fa", repeat_dict=rd, context=context)
    assert context == {
        "trim_start": 501,
        "trim_end": 681,
        "left_trim_method": "exact_anchor",
        "right_trim_method": "exact_anchor",
    }


def test_default_consensus_explicitly_uses_genotype_iupac(tmp_path):
    with patch("muc_one_span.consensus.run_tool", return_value=">c\nAM\n") as run:
        build_consensus(Path("ref.fa"), Path("calls.vcf"), tmp_path / "out.fa")
        command = run.call_args.args[0]
        assert command[command.index("-H") + 1] == "I"
        assert command[command.index("-s") + 1] == "SAMPLE"


def test_trim_exposes_failed_anchor_fallback(tmp_path):
    from muc_one_span.config import load_repeat_dictionary

    full = tmp_path / "full.fa"
    full.write_text(">c\n" + "A" * 1200 + "\n")
    context = {}
    trim_flanking(full, 500, tmp_path / "trim.fa", load_repeat_dictionary(), context=context)
    assert context["left_trim_method"] == "fixed_anchor_not_found"
    assert context["right_trim_method"] == "fixed_anchor_not_found"


def test_identical_trimmed_candidates_do_not_inherit_flank_phase(tmp_path):
    fake = ">c\nAACCGGTT\n"
    alleles = {
        key: {"length": 1, "contig_name": "c", "independent_haplotype_evidence": True}
        for key in ("allele_1", "allele_2")
    }
    with patch("muc_one_span.consensus.run_tool", return_value=fake):
        build_consensus_per_allele(
            tmp_path / "ref.fa",
            dict.fromkeys(alleles, tmp_path / "v.vcf.gz"),
            alleles,
            tmp_path / "out",
            flank_length=2,
        )
    assert all(not a["independent_haplotype_evidence"] for a in alleles.values())
    assert all(a["vntr_phase_status"] == "no_sequence_distinction" for a in alleles.values())
