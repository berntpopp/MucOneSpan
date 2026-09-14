"""Behavior contracts for the maintained batch and timing scripts."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.evaluation.models import (
    Event,
    PredictedAllele,
    RunObservation,
    TruthHaplotype,
    TruthSample,
)


def script_module(name: str):
    path = Path(__file__).parents[2] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"test_{name}_script", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_missing_roots_return_nonzero(tmp_path: Path) -> None:
    missing = tmp_path / "missing"
    assert script_module("benchmark").main(["--data-dir", str(missing)]) != 0
    assert script_module("batch_analyze").main([str(missing), str(tmp_path / "out")]) != 0


def test_expected_sample_without_input_stays_in_denominator(tmp_path: Path) -> None:
    from muc_one_span import benchmarking

    data = tmp_path / "data"
    (data / "missing").mkdir(parents=True)
    truth = TruthSample(
        "missing",
        (TruthHaplotype("h1", "A", ("X",)), TruthHaplotype("h2", "C", ("X",))),
    )
    entries = [{"sample": "missing", "input": "reads.fastq", "platform": "hifi"}]
    records = benchmarking.run_inventory(entries, data, tmp_path / "out")
    with patch("muc_one_span.benchmarking.load_truth", return_value=truth):
        report, failed = benchmarking.evaluate_inventory(entries, data, tmp_path / "out", records)
    assert failed
    assert report["totals"]["truth_haplotypes"] == 2
    assert report["samples"][0]["status"] == "not_attempted"


def test_strict_batch_metrics_reject_wrong_site_and_false_extras(tmp_path: Path) -> None:
    from muc_one_span import benchmarking

    expected = Event(25, "X", "dupC")
    truth = TruthSample("s", (TruthHaplotype("h", "A", ("X",), (expected,)),))
    calls = (
        Event(999, "X", "dupC", True, True, True),
        Event(25, "X", "dupCCCC", True, True, True),
    )
    observation = RunObservation(
        "completed",
        (
            PredictedAllele(
                "p",
                "A",
                ("X",),
                1,
                -8,
                calls,
                sequence_source="p",
                independent_haplotype_evidence=True,
            ),
        ),
        {"exit_code": 0},
    )
    with (
        patch("muc_one_span.benchmarking.load_truth", return_value=truth),
        patch("muc_one_span.benchmarking.load_observation", return_value=observation),
    ):
        report, failed = benchmarking.evaluate_inventory(
            [{"sample": "s"}], tmp_path, tmp_path, [{"sample": "s", "exit_code": 0}]
        )
    assert not failed
    assert report["inventory_mode"] == "explicit"
    assert report["totals"]["metrics"]["event_tp"] == {"min": 0, "max": 0}
    assert report["totals"]["metrics"]["event_fn"] == {"min": 1, "max": 1}
    assert report["totals"]["metrics"]["event_fp"] == {"min": 2, "max": 2}


def test_nonzero_caller_exit_makes_batch_fail(tmp_path: Path) -> None:
    from muc_one_span import benchmarking

    truth = TruthSample("s", (TruthHaplotype("h", "A", ("X",)),))
    observation = RunObservation("insufficient_evidence", run_record={"exit_code": 1})
    record = {"sample": "s", "status": "insufficient_evidence", "exit_code": 1}
    with (
        patch("muc_one_span.benchmarking.load_truth", return_value=truth),
        patch("muc_one_span.benchmarking.load_observation", return_value=observation),
    ):
        _, failed = benchmarking.evaluate_inventory([{"sample": "s"}], tmp_path, tmp_path, [record])
    assert failed


def test_failed_execution_never_scores_stale_output(tmp_path: Path) -> None:
    from muc_one_span.benchmarking import run_pipeline
    from muc_one_span.evaluation import load_observation

    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    output = tmp_path / "output"
    output.mkdir()
    (output / "summary.json").write_text('{"stale": true}')
    with patch("muc_one_span.tools.check_tools", side_effect=RuntimeError("caller failed")):
        record = run_pipeline("sample", reads, output, "hifi", "model", 3)
    assert record["status"] == "execution_failed"
    assert record["exit_code"] != 0
    assert load_observation(output, record).predictions == ()


def test_nonzero_result_overrides_stale_completed_sidecar(tmp_path: Path) -> None:
    from muc_one_span.benchmarking import run_pipeline

    class Result:
        exit_code = 9
        output = "failed"
        exception = RuntimeError("failed")

    class Runner:
        def invoke(self, command, args):
            return Result()

    reads = tmp_path / "reads.fastq"
    reads.touch()
    output = tmp_path / "out"
    output.mkdir()
    (output / "summary.json").write_text('{"stale": true}')
    (output / "run_status.json").write_text('{"schema_version":1,"status":"completed"}')
    record = run_pipeline("sample", reads, output, "hifi", "model", 4, runner=Runner())
    assert record["status"] == "execution_failed"


def test_pipeline_controls_are_forwarded_and_recorded(tmp_path: Path) -> None:
    from muc_one_span.benchmarking import run_pipeline

    class Result:
        exit_code = 0
        output = "ok"
        exception = None

    class Runner:
        def invoke(self, command, args):
            return Result()

    reads = tmp_path / "reads.fq.gz"
    reads.touch()
    model = tmp_path / "models" / "ont"
    record = run_pipeline("sample", reads, tmp_path / "out", "ont", str(model), 7, runner=Runner())
    assert record["platform"] == "ont"
    assert record["model"] == str(model)
    assert record["threads"] == 7
    assert record["cli_args"] == [
        "run",
        "--input",
        str(reads.resolve()),
        "--output-dir",
        str((tmp_path / "out").resolve()),
        "--clair3-model",
        str(model),
        "--threads",
        "7",
        "--platform",
        "ont",
    ]


def test_ambiguous_input_is_an_explicit_failure(tmp_path: Path) -> None:
    from muc_one_span.benchmarking import run_inventory

    sample = tmp_path / "sample"
    sample.mkdir()
    (sample / "reads.bam").touch()
    (sample / "reads.fastq.gz").touch()
    records = run_inventory([{"sample": "sample"}], tmp_path, tmp_path / "out")
    assert records[0]["status"] == "not_attempted"
    assert "ambiguous" in records[0]["error"]


def test_declared_platform_conflict_prevents_caller_run(tmp_path: Path) -> None:
    from muc_one_span.benchmarking import run_inventory

    sample = tmp_path / "sample"
    sample.mkdir()
    (sample / "reads.bam").touch()
    records = run_inventory(
        [{"sample": "sample", "input": "reads.bam", "platform": "ont"}],
        tmp_path,
        tmp_path / "out",
        platform="hifi",
    )
    assert records[0]["status"] == "not_attempted"
    assert "platform conflicts" in records[0]["error"]


def test_simulator_metadata_prevents_silent_platform_mismatch(tmp_path: Path) -> None:
    from muc_one_span.benchmarking import platform_for_sample

    (tmp_path / "reads_metadata.tsv").write_text(
        "Parameter\tValue\nRead_simulation_technology\tONT\n"
    )
    assert platform_for_sample(tmp_path, {}, None) == "ont"
    with pytest.raises(ValueError, match="conflicts with metadata"):
        platform_for_sample(tmp_path, {}, "hifi")


def test_batch_output_declares_legacy_schema_migration(tmp_path: Path) -> None:
    from muc_one_span import benchmarking

    data = tmp_path / "data"
    (data / "sample").mkdir(parents=True)
    inventory = tmp_path / "samples.json"
    inventory.write_text('[{"sample":"sample","input":"missing.bam"}]')
    truth = TruthSample("sample", (TruthHaplotype("h", "A", ("X",)),))
    script = script_module("batch_analyze")
    with patch.object(benchmarking, "load_truth", return_value=truth):
        code = script.main([str(data), str(tmp_path / "out"), "--expected-samples", str(inventory)])
    report = json.loads((tmp_path / "out" / "batch_results.json").read_text())
    assert code != 0
    assert report["schema_version"] == 1
    assert report["schema_migration"]["removed_status"] == "TP_partial"
    assert report["totals"]["truth_haplotypes"] == 1


def test_mixed_platform_inventory_rejects_one_shared_model(tmp_path):
    from muc_one_span.benchmarking import run_inventory

    entries = []
    for platform in ("hifi", "ont"):
        sample = tmp_path / platform
        sample.mkdir()
        (sample / "reads.fq").write_text("@r\nAC\n+\nII\n")
        entries.append({"sample": platform, "platform": platform, "input": "reads.fq"})
    with patch("muc_one_span.benchmarking.run_pipeline") as pipeline:
        records = run_inventory(entries, tmp_path, tmp_path / "out", model="/models/hifi")
    pipeline.assert_not_called()
    assert len(records) == 2
    assert all(r["status"] == "not_attempted" for r in records)
    assert all("per-sample models" in r["error"] for r in records)
