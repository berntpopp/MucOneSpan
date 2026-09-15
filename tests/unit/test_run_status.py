"""Scientific no-call and failed execution must not be indistinguishable."""

import json
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main


def test_pipeline_failure_writes_status_even_when_summary_is_stale(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "summary.json").write_text('{"stale": true}')
    with patch("muc_one_span.tools.check_tools", side_effect=RuntimeError("tool failed")):
        result = CliRunner().invoke(main, ["run", "-i", str(reads), "-o", str(tmp_path)])
    assert result.exit_code != 0
    status = json.loads((tmp_path / "run_status.json").read_text())
    assert status["status"] == "execution_failed"
    assert status["error_type"] == "RuntimeError"


def test_insufficient_evidence_is_typed_but_preserves_valueerror_contract(tmp_path: Path) -> None:
    from muc_one_span.alleles import detect_alleles
    from muc_one_span.run_status import InsufficientEvidenceError, record_run_status

    @record_run_status
    def run(output_dir: str) -> None:
        detect_alleles({51: 2}, min_coverage=10)

    with pytest.raises(InsufficientEvidenceError) as caught:
        run(str(tmp_path))
    assert isinstance(caught.value, ValueError)
    assert (
        json.loads((tmp_path / "run_status.json").read_text())["status"] == "insufficient_evidence"
    )


def test_success_replaces_running_status(tmp_path: Path) -> None:
    from muc_one_span.run_status import record_run_status

    @record_run_status
    def run(output_dir: str) -> int:
        assert json.loads((Path(output_dir) / "run_status.json").read_text())["status"] == "running"
        return 3

    assert run(str(tmp_path)) == 3
    assert json.loads((tmp_path / "run_status.json").read_text())["status"] == "completed"


def test_missing_input_overwrites_stale_completed_status(tmp_path):
    (tmp_path / "run_status.json").write_text('{"schema_version": 1, "status": "completed"}')
    (tmp_path / "summary.json").write_text('{"stale": true}')
    result = CliRunner().invoke(
        main, ["run", "-i", str(tmp_path / "missing.fq"), "-o", str(tmp_path)]
    )
    assert result.exit_code == 2
    status = json.loads((tmp_path / "run_status.json").read_text())
    assert status["status"] == "execution_failed"
    assert status["error_type"] == "BadParameter"


def test_keyboard_interrupt_writes_interrupted_status(tmp_path: Path) -> None:
    from muc_one_span.run_status import record_run_status

    @record_run_status
    def run(output_dir: str) -> None:
        raise KeyboardInterrupt("Simulated user interrupt")

    with pytest.raises(KeyboardInterrupt):
        run(str(tmp_path))

    status = json.loads((tmp_path / "run_status.json").read_text())
    assert status["status"] == "interrupted"
    assert status["error_type"] == "KeyboardInterrupt"


@pytest.mark.parametrize(
    "error,status",
    [(RuntimeError("render failed"), "execution_failed"), (KeyboardInterrupt(), "interrupted")],
)
def test_failed_report_updates_summary_provenance_for_portable_reports(tmp_path, error, status):
    from muc_one_span.run_status import record_run_status

    summary = tmp_path / "summary.json"

    @record_run_status
    def run(output_dir: str) -> None:
        summary.write_text(
            json.dumps({"classifications": {}, "run_status": {"status": "analysis_completed"}})
        )
        raise error

    with pytest.raises(type(error)):
        run(str(tmp_path))
    assert json.loads(summary.read_text())["run_status"]["status"] == status
    assert json.loads((tmp_path / "run_status.json").read_text())["status"] == status


@pytest.mark.parametrize("contents", ["{invalid", "[]"])
def test_malformed_old_summary_does_not_mask_execution_error(tmp_path, contents):
    from muc_one_span.run_status import record_run_status

    (tmp_path / "summary.json").write_text(contents)

    @record_run_status
    def run(output_dir: str) -> None:
        raise RuntimeError("original execution error")

    with pytest.raises(RuntimeError, match="original execution error"):
        run(str(tmp_path))
    assert json.loads((tmp_path / "run_status.json").read_text())["status"] == "execution_failed"
