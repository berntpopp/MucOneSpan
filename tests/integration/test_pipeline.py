"""End-to-end pipeline tests using MucOneUp-generated test data."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from tests.conftest import (
    requires_bcftools,
    requires_clair3,
    requires_minimap2,
    requires_samtools,
)

TESTDATA_DIR = Path(__file__).resolve().parents[1] / "data" / "generated"


def sample_input(sample: str) -> Path:
    """Resolve a sample independently, so partially generated datasets can run."""
    sample_dir = TESTDATA_DIR / sample
    for pattern in ("*.bam", "*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"):
        inputs = sorted(sample_dir.glob(pattern))
        if inputs:
            return inputs[0]
    pytest.skip(f"MucOneUp reads unavailable for {sample}; run 'make generate-testdata'")


@requires_minimap2
@requires_samtools
@pytest.mark.e2e
class TestAlleleLengthDetection:
    """Test allele length detection against ground truth."""

    @pytest.mark.parametrize(
        "sample,expected_h1,expected_h2",
        [
            ("sample_dupc_60_80", 60, 80),
            ("sample_normal_60_80", 60, 80),
            pytest.param(
                "sample_homozygous_60_60",
                60,
                60,
                marks=pytest.mark.xfail(
                    strict=True,
                    raises=AssertionError,
                    reason="Known limitation: indel-valley splitting separates same-length alleles (docs/reference/limitations.md)",
                ),
            ),
            pytest.param(
                "sample_asymmetric_25_140",
                25,
                140,
                marks=pytest.mark.xfail(
                    strict=True,
                    raises=AssertionError,
                    reason="Known limitation: extreme PCR bias obscures the 140-repeat allele (docs/reference/limitations.md)",
                ),
            ),
            ("sample_short_25_30", 25, 30),
            ("sample_long_120_140", 120, 140),
        ],
    )
    def test_allele_lengths(
        self, sample: str, expected_h1: int, expected_h2: int, tmp_path: Path
    ) -> None:
        """Detected allele lengths match ground truth within tolerance."""
        from muc_one_span.alleles import detect_alleles, parse_idxstats
        from muc_one_span.config import load_repeat_dictionary
        from muc_one_span.ladder import generate_ladder_fasta
        from muc_one_span.mapping import get_idxstats, map_reads

        input_path = sample_input(sample)
        rd = load_repeat_dictionary()
        ref = tmp_path / "ladder.fa"
        generate_ladder_fasta(rd, ref)

        bam = map_reads(input_path, ref, tmp_path)
        idxstats = get_idxstats(bam)
        counts = parse_idxstats(idxstats)
        result = detect_alleles(counts, min_coverage=10, bam_path=bam)

        detected = sorted([result["allele_1"]["length"], result["allele_2"]["length"]])
        expected = sorted([expected_h1, expected_h2])

        # Allow +/- 2 repeat tolerance
        assert abs(detected[0] - expected[0]) <= 2, f"Expected {expected[0]}, got {detected[0]}"
        assert abs(detected[1] - expected[1]) <= 2, f"Expected {expected[1]}, got {detected[1]}"


@requires_minimap2
@requires_samtools
@requires_bcftools
@requires_clair3
@pytest.mark.e2e
class TestFullPipeline:
    """Full pipeline e2e tests with mutation detection validation."""

    @pytest.mark.parametrize(
        "sample,expect_mutation",
        [
            ("sample_dupc_60_80", True),
            ("sample_normal_60_80", False),
        ],
    )
    def test_full_pipeline(
        self, tmp_path: Path, sample: str, expect_mutation: bool, clair3_model: Path
    ) -> None:
        """Execute every CLI stage and check VCF-backed frameshift calls."""
        from muc_one_span.cli import main
        from muc_one_span.config import load_repeat_dictionary
        from muc_one_span.ladder import generate_ladder_fasta

        input_path = sample_input(sample)
        reference = tmp_path / "ladder.fa"
        generate_ladder_fasta(load_repeat_dictionary(), reference)
        output = tmp_path / "pipeline"
        result = CliRunner().invoke(
            main,
            [
                "run",
                "--input",
                str(input_path),
                "--reference",
                str(reference),
                "--output-dir",
                str(output),
                "--clair3-model",
                str(clair3_model),
                "--threads",
                "2",
            ],
        )
        assert result.exit_code == 0, f"{result.output}\n{result.exception}"
        summary = json.loads((output / "summary.json").read_text())
        classifications = summary["classifications"]
        assert classifications
        assert all(value["structure"] for value in classifications.values())
        frameshifts = [
            mutation
            for value in classifications.values()
            for mutation in value["mutations"]
            if mutation.get("frameshift") and mutation.get("vcf_support")
        ]
        if expect_mutation:
            assert frameshifts, summary
            assert any(mutation.get("mutation_name") == "dupC" for mutation in frameshifts)
        else:
            assert not frameshifts, summary
