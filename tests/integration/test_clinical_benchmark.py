"""Bounded real-process validation of the clinical worker transport boundary."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from muc_one_span.tools import run_tool_pipeline


@pytest.mark.integration
def test_clinical_worker_preserves_cli_failure_separately_from_transport(tmp_path: Path) -> None:
    output = tmp_path / "attempt"
    output.mkdir()
    config = tmp_path / "invocation.json"
    config.write_text(
        json.dumps(
            {
                "output": str(output),
                "argv": [
                    "run",
                    "--input",
                    str(tmp_path / "missing.fastq"),
                    "--output-dir",
                    str(output),
                    "--clair3-model",
                    str(tmp_path / "missing-model"),
                ],
            }
        )
    )
    run_tool_pipeline(
        [[sys.executable, "-m", "muc_one_span.clinical_worker", str(config)]], timeout=10
    )
    worker = json.loads((output / "worker.json").read_text())
    assert worker["exit_code"] == 2
    assert "does not exist" in worker["error"].lower()
    assert (output / "cli.log").is_file()
    assert not (output / "summary.json").exists()
