"""The ladder engine is deprecated: a run warns and summary.json records it (additive)."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_span.cli import main
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.deprecations import LADDER_DEPRECATION_MESSAGE, run_deprecations
from muc_one_span.pipeline_tail import finish_run
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings

_CLS = {"structure": "1 2 3 X 6 7 8 9", "repeats": [], "mutations_detected": []}


def _with_engine(engine: str) -> RuntimeSettings:
    return dataclasses.replace(
        DEFAULT_SETTINGS, run=dataclasses.replace(DEFAULT_SETTINGS.run, engine=engine)
    )


def test_default_engine_has_no_deprecation() -> None:
    assert run_deprecations(DEFAULT_SETTINGS) == []


def test_ladder_engine_is_recorded_as_deprecated() -> None:
    (record,) = run_deprecations(_with_engine("ladder"))
    assert (record["setting"], record["value"], record["status"]) == (
        "run.engine",
        "ladder",
        "deprecated",
    )
    assert record["replacement"] == DEFAULT_SETTINGS.run.engine
    assert record["message"] == LADDER_DEPRECATION_MESSAGE
    assert "--engine ladder" in record["message"] and "future release" in record["message"]


def _finish(tmp_path: Path, settings: RuntimeSettings) -> dict:
    fa = tmp_path / "a1.fa"
    fa.write_text(">a1\nACGT\n")
    with patch("muc_one_span.classify.classify_sequence", return_value=_CLS):
        finish_run(
            out=tmp_path,
            input_path=str(tmp_path / "reads.fastq"),
            rd=load_repeat_dictionary(None),
            settings=settings,
            alleles_result={"allele_1": {"length": 1}},
            consensus_paths={"allele_1": fa},
            vcf_paths={},
            tool_versions={},
            configuration_record={},
            report=False,
            bam_path=None,
            fasta_path=None,
        )
    return json.loads((tmp_path / "summary.json").read_text())


def test_summary_records_deprecations_for_both_engines(tmp_path: Path) -> None:
    (tmp_path / "h").mkdir()
    (tmp_path / "l").mkdir()
    assert _finish(tmp_path / "h", DEFAULT_SETTINGS)["deprecations"] == []
    ladder = _finish(tmp_path / "l", _with_engine("ladder"))["deprecations"]
    assert [d["value"] for d in ladder] == ["ladder"]


def test_ladder_run_warns_on_stderr_and_records_it(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.touch()
    out = tmp_path / "results"
    fa = tmp_path / "a1.fa"
    fa.write_text(">a1\nACGT\n")
    alleles = {"homozygous": True, "allele_1": {"length": 60, "contig_name": "c51"}}
    with (
        patch("muc_one_span.tools.check_tools"),
        patch("muc_one_span.tools.get_tool_versions", return_value={}),
        patch("muc_one_span.mapping.map_reads", return_value=tmp_path / "m.bam"),
        patch("muc_one_span.mapping.get_idxstats", return_value="c51\t3060\t100\t0\n"),
        patch("muc_one_span.alleles.parse_idxstats", return_value={51: 100}),
        patch("muc_one_span.alleles.detect_alleles", return_value=alleles),
        patch("muc_one_span.calling.call_variants_per_allele", return_value={}),
        patch("muc_one_span.consensus.build_consensus_per_allele", return_value={"allele_1": fa}),
        patch("muc_one_span.classify.classify_sequence", return_value=_CLS),
    ):
        result = CliRunner().invoke(
            main, ["run", "--input", str(reads), "--output-dir", str(out), "--engine", "ladder"]
        )
    assert result.exit_code == 0, result.output
    assert LADDER_DEPRECATION_MESSAGE in " ".join(result.stderr.split())
    summary = json.loads((out / "summary.json").read_text())
    assert [d["value"] for d in summary["deprecations"]] == ["ladder"]
