"""``run_split`` over a manifest, with the real Task 6 manifest field names."""

from __future__ import annotations

import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.benchsim.generate import FASTQ
from muc_one_span.benchsim.run_cases import platform_for_profile, run_split


def _row(design_id: str, profile: str, status: str, **extra: object) -> dict[str, object]:
    """A manifest row shaped like ``generate.generate_case`` actually writes."""
    return {
        "design_id": design_id,
        "split": "dev",
        "profile": profile,
        "design": {"profile": profile},
        "muconeup_version": "0.45.0",
        "status": status,
        "error": extra.pop("error", None),
        **extra,
    }


def _manifest(split_dir: Path, rows: list[dict[str, object]]) -> Path:
    split_dir.mkdir(parents=True, exist_ok=True)
    path = split_dir / "manifest.jsonl"
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))
    return path


def _write_reads(split_dir: Path, design_id: str, profile: str) -> None:
    reads_dir = split_dir / design_id / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    fastq = reads_dir / FASTQ[profile].format(design_id)
    fastq.write_text("@r\nACGT\n+\nIIII\n")


def _write_truth(split_dir: Path, design_id: str) -> None:
    (split_dir / design_id / "truth").mkdir(parents=True, exist_ok=True)


def _r1041(platform: str) -> str:
    """Module-level (picklable) ``model_for`` stand-in, for the process-pool test."""
    return "r1041"


def test_denominator_kept_and_engine_forwarded(tmp_path: Path) -> None:
    """An ``ok`` row runs on every engine; a non-``ok`` row keeps the denominator."""
    split_dir = tmp_path / "dev"
    _write_reads(split_dir, "dev-a", "ont_amplicon_r10")
    _write_truth(split_dir, "dev-a")
    rows = [
        _row("dev-a", "ont_amplicon_r10", "ok"),
        _row("dev-b", "hifi_amplicon", "design_invalid", error="truth events != design"),
    ]
    manifest = _manifest(split_dir, rows)

    calls = []

    def fake(
        sample, input_path, output_dir, platform, model, threads, *, runner=None, engine="ladder"
    ):
        calls.append((sample, platform, engine))
        return {"sample": sample, "status": "completed", "result_dir": str(output_dir)}

    with patch("muc_one_span.benchsim.run_cases.run_pipeline", side_effect=fake):
        records = run_split(
            manifest, ["ladder", "hybrid"], tmp_path / "res", lambda p: "r1041", 1, 1
        )

    assert sorted(calls) == [("dev-a", "ont", "hybrid"), ("dev-a", "ont", "ladder")]
    statuses = {(r["engine"], r["sample"]): r["status"] for r in records}
    assert statuses[("ladder", "dev-a")] == "completed"
    assert statuses[("hybrid", "dev-a")] == "completed"
    assert statuses[("ladder", "dev-b")] == "not_attempted"
    assert statuses[("hybrid", "dev-b")] == "not_attempted"
    assert len(records) == 4

    ladder_measurements = json.loads((tmp_path / "res/ladder/measurements.json").read_text())
    assert ladder_measurements
    assert {r["sample"] for r in ladder_measurements} == {"dev-a", "dev-b"}


def test_not_attempted_row_never_reaches_the_caller(tmp_path: Path) -> None:
    """``design_invalid``/``generation_failed`` rows must not call ``run_pipeline``."""
    split_dir = tmp_path / "dev"
    rows = [_row("dev-b", "hifi_amplicon", "generation_failed", error="MucOneUp crashed")]
    manifest = _manifest(split_dir, rows)

    with patch("muc_one_span.benchsim.run_cases.run_pipeline") as pipeline:
        records = run_split(manifest, ["ladder"], tmp_path / "res", lambda p: "m", 1, 1)

    pipeline.assert_not_called()
    assert records[0]["status"] == "not_attempted"
    assert records[0]["error"] == "MucOneUp crashed"
    assert records[0]["platform"] == "hifi"


def test_run_pipeline_status_is_preserved_when_output_is_invalid(tmp_path: Path) -> None:
    """A run that was attempted but produced no ``summary.json`` keeps its real status."""
    split_dir = tmp_path / "dev"
    _write_reads(split_dir, "dev-c", "ont_amplicon_r10")
    _write_truth(split_dir, "dev-c")
    manifest = _manifest(split_dir, [_row("dev-c", "ont_amplicon_r10", "ok")])

    def fake(
        sample, input_path, output_dir, platform, model, threads, *, runner=None, engine="ladder"
    ):
        return {
            "sample": sample,
            "status": "invalid_artifacts",
            "result_dir": str(output_dir),
            "error": "successful command did not produce summary.json",
        }

    with patch("muc_one_span.benchsim.run_cases.run_pipeline", side_effect=fake):
        records = run_split(manifest, ["ladder"], tmp_path / "res", lambda p: "r1041", 1, 1)

    assert records[0]["status"] == "invalid_artifacts"
    assert records[0]["engine"] == "ladder"


def test_inventory_has_truth_and_result_dir_per_case(tmp_path: Path) -> None:
    split_dir = tmp_path / "dev"
    _write_reads(split_dir, "dev-a", "ont_amplicon_r10")
    _write_truth(split_dir, "dev-a")
    manifest = _manifest(split_dir, [_row("dev-a", "ont_amplicon_r10", "ok")])

    def fake(
        sample, input_path, output_dir, platform, model, threads, *, runner=None, engine="ladder"
    ):
        return {"sample": sample, "status": "completed", "result_dir": str(output_dir)}

    with patch("muc_one_span.benchsim.run_cases.run_pipeline", side_effect=fake):
        run_split(manifest, ["ladder"], tmp_path / "res", lambda p: "r1041", 1, 1)

    inventory = json.loads((tmp_path / "res/ladder/inventory_dev.json").read_text())
    assert inventory == [
        {
            "sample": "dev-a",
            "truth_dir": str(split_dir / "dev-a" / "truth"),
            "result_dir": str(tmp_path / "res" / "ladder" / "dev-a"),
        }
    ]


def test_platform_for_profile_mapping() -> None:
    assert platform_for_profile("ont_amplicon_r10") == "ont"
    assert platform_for_profile("ont_genomic_targeted") == "ont"
    assert platform_for_profile("hifi_amplicon") == "hifi"
    with pytest.raises(ValueError, match="no platform mapping"):
        platform_for_profile("unknown_profile")
    with pytest.raises(TypeError):
        platform_for_profile(None)


def test_unmapped_profile_becomes_not_attempted_without_calling_the_caller(tmp_path: Path) -> None:
    split_dir = tmp_path / "dev"
    manifest = _manifest(split_dir, [_row("dev-z", "bogus_profile", "ok")])

    with patch("muc_one_span.benchsim.run_cases.run_pipeline") as pipeline:
        records = run_split(manifest, ["ladder"], tmp_path / "res", lambda p: "m", 1, 1)

    pipeline.assert_not_called()
    assert records[0]["status"] == "not_attempted"
    assert "no platform mapping" in records[0]["error"]


def test_run_pipeline_exception_becomes_execution_failed(tmp_path: Path) -> None:
    """A caller crash (e.g. a missing external tool) keeps the denominator, never propagates."""
    split_dir = tmp_path / "dev"
    _write_reads(split_dir, "dev-a", "ont_amplicon_r10")
    _write_truth(split_dir, "dev-a")
    manifest = _manifest(split_dir, [_row("dev-a", "ont_amplicon_r10", "ok")])

    with patch(
        "muc_one_span.benchsim.run_cases.run_pipeline", side_effect=RuntimeError("tool missing")
    ):
        records = run_split(manifest, ["ladder"], tmp_path / "res", lambda p: "r1041", 1, 1)

    assert records[0]["status"] == "execution_failed"  # attempted, then crashed
    assert records[0]["error"] == "tool missing"
    assert records[0]["engine"] == "ladder"


def test_ok_row_without_a_reads_file_becomes_not_attempted(tmp_path: Path) -> None:
    """A trusted ``ok`` row whose reads went missing keeps the denominator, never crashes."""
    split_dir = tmp_path / "dev"
    manifest = _manifest(split_dir, [_row("dev-a", "ont_amplicon_r10", "ok")])

    with patch("muc_one_span.benchsim.run_cases.run_pipeline") as pipeline:
        records = run_split(manifest, ["ladder"], tmp_path / "res", lambda p: "r1041", 1, 1)

    pipeline.assert_not_called()
    assert records[0]["status"] == "not_attempted"
    assert "reads file not found" in records[0]["error"]


def test_jobs_parallel_matches_sequential(tmp_path: Path) -> None:
    """The process-pool path (``jobs`` > 1) aggregates the same records as sequential."""
    split_dir = tmp_path / "dev"
    _write_reads(split_dir, "dev-a", "ont_amplicon_r10")
    _write_reads(split_dir, "dev-b", "hifi_amplicon")
    _write_truth(split_dir, "dev-a")
    _write_truth(split_dir, "dev-b")
    rows = [
        _row("dev-a", "ont_amplicon_r10", "ok"),
        _row("dev-b", "hifi_amplicon", "ok"),
    ]
    manifest = _manifest(split_dir, rows)

    def fake(
        sample, input_path, output_dir, platform, model, threads, *, runner=None, engine="ladder"
    ):
        return {"sample": sample, "status": "completed", "result_dir": str(output_dir)}

    # Fork workers so the patched caller reaches them on every Python (3.14 defaults
    # to forkserver, whose workers re-import the real module).
    fork_pool = partial(ProcessPoolExecutor, mp_context=multiprocessing.get_context("fork"))
    with (
        patch("muc_one_span.benchsim.run_cases.run_pipeline", side_effect=fake),
        patch("muc_one_span.benchsim.run_cases.ProcessPoolExecutor", fork_pool),
    ):
        records = run_split(manifest, ["ladder"], tmp_path / "res", _r1041, 1, 2)

    statuses = {r["sample"]: r["status"] for r in records}
    assert statuses == {"dev-a": "completed", "dev-b": "completed"}


def test_config_is_forwarded_only_when_given(tmp_path: Path) -> None:
    split_dir = tmp_path / "dev"
    _write_reads(split_dir, "dev-a", "ont_amplicon_r10")
    _write_truth(split_dir, "dev-a")
    manifest = _manifest(split_dir, [_row("dev-a", "ont_amplicon_r10", "ok")])
    seen: list[dict[str, object]] = []

    def fake(sample, input_path, output_dir, platform, model, threads, **kwargs):
        seen.append(kwargs)
        return {"sample": sample, "status": "completed", "result_dir": str(output_dir)}

    config = tmp_path / "overlay.json"
    with patch("muc_one_span.benchsim.run_cases.run_pipeline", side_effect=fake):
        run_split(manifest, ["ladder"], tmp_path / "a", lambda p: "m", 1, 1)
        run_split(manifest, ["ladder"], tmp_path / "b", lambda p: "m", 1, 1, config=config)
    assert "config" not in seen[0]
    assert seen[1]["config"] == config


def test_model_lookup_failure_is_not_attempted(tmp_path: Path) -> None:
    """No model for the platform: the caller was never invoked."""
    split_dir = tmp_path / "dev"
    _write_reads(split_dir, "dev-a", "ont_amplicon_r10")
    _write_truth(split_dir, "dev-a")
    manifest = _manifest(split_dir, [_row("dev-a", "ont_amplicon_r10", "ok")])

    def no_model(platform: str) -> str:
        raise KeyError(platform)

    with patch("muc_one_span.benchsim.run_cases.run_pipeline") as caller:
        records = run_split(manifest, ["ladder"], tmp_path / "res", no_model, 1, 1)
    caller.assert_not_called()
    assert records[0]["status"] == "not_attempted"
