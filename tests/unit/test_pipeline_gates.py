"""Pipeline wiring of pre-analysis checks and allele-selection gates (all tools mocked)."""

from __future__ import annotations

import json
from collections.abc import Iterable
from contextlib import AbstractContextManager, ExitStack
from pathlib import Path
from typing import Any
from unittest.mock import patch

from click.testing import CliRunner, Result

from muc_one_span.cli import main


def _alleles() -> dict[str, Any]:
    def allele(units: int, primary: int) -> dict[str, Any]:
        return {
            "length": units + 9,
            "reference_length": units + 9,
            "reads": primary * 6,
            "primary_alignment_records": primary,
            "canonical_repeats": units,
            "contig_name": f"contig_{units}",
            "cluster_contigs": [f"contig_{units}"],
            "fit_metrics": {f"contig_{units}": {"primary_alignment_records": primary}},
        }

    return {
        "homozygous": False,
        "same_length": False,
        "allele_1": allele(51, 12),
        "allele_2": allele(71, 90),
    }


def run_mocked_pipeline(
    tmp_path: Path,
    *args: str,
    extra_patches: Iterable[AbstractContextManager[Any]] = (),
) -> tuple[Result, list[str]]:
    """Invoke ``run`` with every tool stage mocked; record check/map call order."""
    input_file = tmp_path / "reads.fastq"
    input_file.touch()
    fasta = tmp_path / "allele.fa"
    fasta.write_text(">allele_vntr\nACGT\n")
    classification = {
        "structure": "1 2 3 X 6 7 8 9",
        "repeats": [],
        "mutations_detected": [],
        "allele_confidence": 1.0,
    }
    calls: list[str] = []

    def check_tools(tools: list[str]) -> bool:
        calls.append("check:" + ",".join(tools))
        return True

    def map_reads(*_args: Any, **_kwargs: Any) -> Path:
        calls.append("map")
        return tmp_path / "mapping.bam"

    with ExitStack() as stack:
        for manager in (
            patch("muc_one_span.tools.check_tools", side_effect=check_tools),
            patch("muc_one_span.tools.get_tool_versions", return_value={}),
            patch("muc_one_span.mapping.map_reads", side_effect=map_reads),
            patch("muc_one_span.mapping.get_idxstats", return_value="contig_51\t3120\t72\t0\n"),
            patch("muc_one_span.alleles.parse_idxstats", return_value={51: 72, 71: 540}),
            patch("muc_one_span.alleles.detect_alleles", side_effect=lambda *a, **k: _alleles()),
            patch(
                "muc_one_span.calling.call_variants_per_allele",
                return_value={
                    "allele_1": tmp_path / "a1.vcf.gz",
                    "allele_2": tmp_path / "a2.vcf.gz",
                },
            ),
            patch(
                "muc_one_span.consensus.build_consensus_per_allele",
                return_value={"allele_1": fasta, "allele_2": fasta},
            ),
            patch("muc_one_span.classify.classify_sequence", return_value=classification),
            patch(
                "muc_one_span.classify.validate_mutations_against_vcf", return_value=classification
            ),
            patch("muc_one_span.vcf.parse_vcf_variants", return_value=[]),
            *extra_patches,
        ):
            stack.enter_context(manager)
        result = CliRunner().invoke(
            main, ["run", "--input", str(input_file), "--output-dir", str(tmp_path / "out"), *args]
        )
    return result, calls


def test_alleles_json_records_selection_and_depth_gates(tmp_path: Path) -> None:
    result, _ = run_mocked_pipeline(tmp_path)
    assert result.exit_code == 0, result.output
    alleles = json.loads((tmp_path / "out" / "alleles.json").read_text())
    assert alleles["allele_1"]["depth_status"] == "low"
    assert alleles["allele_1"]["depth_threshold"] == 30
    assert alleles["allele_2"]["depth_status"] == "adequate"
    assert alleles["allele_2"]["selection_status"] == "resolved"
    assert alleles["allele_2"]["length_status"] == "consistent_with_consensus_contig"


def test_failed_igv_preflight_stops_before_mapping(tmp_path: Path) -> None:
    failing = patch(
        "muc_one_span.report_igv.preflight_igv_report",
        side_effect=RuntimeError("IGV report preflight failed for --report-igv embedded"),
    )
    result, calls = run_mocked_pipeline(
        tmp_path, "--report-igv", "embedded", extra_patches=[failing]
    )
    assert result.exit_code != 0
    assert "map" not in calls
    assert calls == ["check:minimap2,samtools,bcftools,run_clair3.sh,create_report"]
    status = json.loads((tmp_path / "out" / "run_status.json").read_text())
    assert status["status"] == "execution_failed"


def test_igv_off_skips_preflight_and_create_report(tmp_path: Path) -> None:
    probe = patch("muc_one_span.report_igv.preflight_igv_report")
    with probe as preflight:
        result, calls = run_mocked_pipeline(tmp_path)
    assert result.exit_code == 0, result.output
    preflight.assert_not_called()
    assert calls[0] == "check:minimap2,samtools,bcftools,run_clair3.sh"
