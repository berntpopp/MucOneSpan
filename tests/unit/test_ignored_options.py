"""Ladder-only options given to a hybrid run warn, and are recorded as ignored (additive)."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main
from muc_one_span.deprecations import LADDER_ONLY_OPTIONS
from muc_one_span.pipeline import execute_pipeline
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings

WARNING = "is ignored by the hybrid engine; use --engine ladder (deprecated)"


def _invoke(tmp_path: Path, *extra: str, config: dict | None = None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    prefix: list[str] = []
    if config is not None:
        path = tmp_path / "settings.json"
        path.write_text(json.dumps({"schema_version": 1, **config}))
        prefix = ["--config", str(path)]
    args = [*prefix, "run", "-i", str(reads), "-o", str(tmp_path / "out"), *extra]
    with patch("muc_one_span.pipeline._run_hybrid") as hybrid:
        result = CliRunner().invoke(main, args)
    return result, hybrid


def _record(tmp_path: Path) -> dict:
    return json.loads((tmp_path / "out" / "run_configuration.json").read_text())


HYBRID_UNUSED_FLAGS = {
    "--clair3-model",
    "--min-qual",
    "--minimap2-preset",
    "--platform",
    "--min-coverage",
    "--threads",
    "--mapping-timeout",
    "--reference",
}


def test_ladder_only_option_map_names_every_option_the_hybrid_path_ignores() -> None:
    assert set(LADDER_ONLY_OPTIONS.values()) == HYBRID_UNUSED_FLAGS


def test_every_run_option_is_either_used_by_hybrid_or_listed() -> None:
    from muc_one_span.cli_run import run

    flags = {opt for p in run.params for opt in p.opts if opt.startswith("--")}
    used = {"--input", "--output-dir", "--report", "--report-igv", "--engine", "--assay"}
    used |= {"--igv-session"}
    assert flags - used == HYBRID_UNUSED_FLAGS


def test_common_ladder_options_warn_on_a_hybrid_run(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">c\nACGT\n")
    args = ["--platform", "ont", "--min-coverage", "50", "--threads", "8"]
    args += ["--mapping-timeout", "60", "--reference", str(ref)]
    result, hybrid = _invoke(tmp_path, *args)
    assert result.exit_code == 0, result.output
    stderr = " ".join(result.stderr.split())
    for flag in ("--platform", "--min-coverage", "--threads", "--mapping-timeout", "--reference"):
        assert f"Warning: {flag} {WARNING}" in stderr
    hybrid.assert_called_once()


def test_custom_dictionary_does_not_require_a_reference_for_hybrid(tmp_path: Path) -> None:
    layout = {"pre": ["1", "2", "3", "4"], "after": ["6", "7", "8", "9"]}
    result, hybrid = _invoke(tmp_path, config={"reference_layout": layout})
    assert result.exit_code == 0, result.output
    hybrid.assert_called_once()
    ladder, _ = _invoke(tmp_path / "l", "--engine", "ladder", config={"reference_layout": layout})
    assert ladder.exit_code == 2 and "--reference" in ladder.output


def test_explicit_ladder_options_warn_and_are_recorded_as_ignored(tmp_path: Path) -> None:
    model = tmp_path / "model"
    model.mkdir()
    result, hybrid = _invoke(
        tmp_path, "--clair3-model", str(model), "--min-qual", "5", "--minimap2-preset", "map-ont"
    )
    assert result.exit_code == 0, result.output
    stderr = " ".join(result.stderr.split())
    for flag in ("--clair3-model", "--min-qual", "--minimap2-preset"):
        assert f"Warning: {flag} {WARNING}" in stderr
    record = _record(tmp_path)
    assert [o["option"] for o in record["ignored_options"]] == [
        "--clair3-model",
        "--min-qual",
        "--minimap2-preset",
    ]
    assert record["model_selection"] == "not used (hybrid engine)"
    assert record["resolved_minimap2_preset"] is None
    assert hybrid.call_args.args[4]["ignored_options"] == record["ignored_options"]


def test_default_hybrid_run_ignores_nothing(tmp_path: Path) -> None:
    result, _ = _invoke(tmp_path)
    assert result.exit_code == 0, result.output
    assert "ignored by the hybrid engine" not in result.stderr
    assert _record(tmp_path)["ignored_options"] == []


def test_non_default_config_value_is_ignored_but_config_defaults_are_not(tmp_path: Path) -> None:
    # A config file re-applies every run value through Click's default_map; only the
    # non-default ladder-only value counts, and hybrid-used values (assay) never do.
    result, _ = _invoke(tmp_path, config={"run": {"min_qual": 20, "assay": "genomic"}})
    assert result.exit_code == 0, result.output
    assert f"Warning: --min-qual {WARNING}" in " ".join(result.stderr.split())
    ignored = _record(tmp_path)["ignored_options"]
    assert [(o["setting"], o["value"]) for o in ignored] == [("run.min_qual", 20)]


def test_ladder_run_ignores_nothing_and_keeps_model_provenance(tmp_path: Path) -> None:
    model = tmp_path / "model"
    model.mkdir()
    with patch("muc_one_span.tools.check_tools", side_effect=RuntimeError("stop after config")):
        result, _ = _invoke(tmp_path, "--engine", "ladder", "--clair3-model", str(model))
    assert "ignored by the hybrid engine" not in result.stderr
    record = _record(tmp_path)
    assert record["ignored_options"] == []
    assert record["model_selection"] == "explicit path"
    assert record["resolved_minimap2_preset"] == "map-hifi"


def test_report_igv_with_hybrid_runs_and_uses_threads_and_mapping_timeout(
    tmp_path: Path,
) -> None:
    result, hybrid = _invoke(
        tmp_path, "--report-igv", "embedded", "--threads", "2", "--mapping-timeout", "60"
    )
    assert result.exit_code == 0, result.output
    hybrid.assert_called_once()
    assert "ignored by the hybrid engine" not in result.stderr
    assert _record(tmp_path)["ignored_options"] == []


def test_hybrid_without_igv_still_ignores_threads(tmp_path: Path) -> None:
    result, _ = _invoke(tmp_path, "--threads", "2")
    assert result.exit_code == 0, result.output
    assert [r["option"] for r in _record(tmp_path)["ignored_options"]] == ["--threads"]


def _execute(tmp_path: Path, settings: RuntimeSettings, **kwargs: str) -> None:
    """``execute_pipeline`` with no ``engine`` argument: the engine comes from settings."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    args = (str(reads), str(tmp_path / "out"), None, "", 1, 10, 5.0, False, "ont", None)
    execute_pipeline(*args, settings=settings, **kwargs)


def test_engine_none_takes_the_engine_from_the_configuration(tmp_path: Path) -> None:
    hybrid = replace(DEFAULT_SETTINGS, run=replace(DEFAULT_SETTINGS.run, engine="hybrid"))
    with patch("muc_one_span.pipeline._run_hybrid") as run_hybrid:
        _execute(tmp_path / "h", hybrid)
    assert run_hybrid.call_count == 1

    class LadderPathError(Exception):
        pass

    ladder = replace(DEFAULT_SETTINGS, run=replace(DEFAULT_SETTINGS.run, engine="ladder"))
    with (
        patch("muc_one_span.pipeline._run_hybrid") as run_hybrid,
        patch("muc_one_span.tools.check_tools", side_effect=LadderPathError),
        pytest.raises(LadderPathError),
    ):
        _execute(tmp_path / "l", ladder)
    assert run_hybrid.call_count == 0


def test_config_only_hybrid_engine_with_igv_reaches_the_hybrid_path(tmp_path: Path) -> None:
    run = replace(DEFAULT_SETTINGS.run, engine="hybrid", report_igv="embedded")
    configured = replace(DEFAULT_SETTINGS, run=run)
    with patch("muc_one_span.pipeline._run_hybrid") as run_hybrid:
        _execute(tmp_path, configured, report_igv="embedded")
    assert run_hybrid.call_count == 1


def test_igv_session_uses_threads_and_is_refused_by_the_ladder(tmp_path: Path) -> None:
    result, hybrid = _invoke(tmp_path / "h", "--igv-session", "--threads", "2")
    assert result.exit_code == 0, result.output
    hybrid.assert_called_once()
    record = json.loads((tmp_path / "h" / "out" / "run_configuration.json").read_text())
    assert record["ignored_options"] == []
    assert record["settings"]["run"]["igv_session"] is True
    ladder, _ = _invoke(tmp_path / "l", "--engine", "ladder", "--igv-session")
    assert ladder.exit_code == 2 and "hybrid engine only" in " ".join(ladder.output.split())
