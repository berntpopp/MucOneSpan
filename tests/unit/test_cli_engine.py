"""The run command's --engine/--assay options and their configuration defaults."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_span.cli import main


def _fastq(tmp_path: Path) -> Path:
    fq = tmp_path / "r.fastq"
    fq.write_text("@r\nACGT\n+\nIIII\n")
    return fq


def test_engine_option_reaches_pipeline(tmp_path: Path) -> None:
    args = ["run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path / "o")]
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, [*args, "--engine", "hybrid", "--assay", "genomic"])
    assert res.exit_code == 0, res.output
    assert (ex.call_args.kwargs["engine"], ex.call_args.kwargs["assay"]) == ("hybrid", "genomic")


def test_engine_defaults_to_hybrid(tmp_path: Path) -> None:
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, ["run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path)])
    assert res.exit_code == 0, res.output
    assert (ex.call_args.kwargs["engine"], ex.call_args.kwargs["assay"]) == ("hybrid", "amplicon")


def test_ladder_is_still_selectable_from_config(tmp_path: Path) -> None:
    cfg = tmp_path / "c.json"
    cfg.write_text('{"schema_version": 1, "run": {"engine": "ladder"}}')
    args = ["--config", str(cfg), "run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path / "o")]
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, args)
    assert res.exit_code == 0, res.output
    assert ex.call_args.kwargs["engine"] == "ladder"


def test_engine_help_names_hybrid_default_and_ladder_deprecation() -> None:
    res = CliRunner().invoke(main, ["run", "--help"])
    assert res.exit_code == 0, res.output
    text = " ".join(res.output.split())
    assert "default: hybrid" in text and "deprecated" in text


def test_engine_default_comes_from_config(tmp_path: Path) -> None:
    cfg = tmp_path / "c.json"
    cfg.write_text('{"schema_version": 1, "run": {"engine": "hybrid"}}')
    args = ["--config", str(cfg), "run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path / "o")]
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, args)
    assert res.exit_code == 0, res.output
    assert ex.call_args.kwargs["engine"] == "hybrid"


def test_run_is_still_importable_from_cli() -> None:
    from muc_one_span.cli import run
    from muc_one_span.cli_run import run as moved

    assert run is moved and main.commands["run"] is moved
