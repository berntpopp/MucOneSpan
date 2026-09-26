"""Ladder-only options given to a hybrid run warn, and are recorded as ignored (additive)."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_span.cli import main
from muc_one_span.deprecations import LADDER_ONLY_OPTIONS

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


def test_report_igv_with_hybrid_fails_before_configuration_names_both_remedies(
    tmp_path: Path,
) -> None:
    result, hybrid = _invoke(tmp_path, "--report-igv", "embedded")
    assert result.exit_code == 2
    text = " ".join(result.output.split())
    assert "--engine ladder (deprecated)" in text and "--report-igv off" in text
    assert not (tmp_path / "out" / "run_configuration.json").exists()
    hybrid.assert_not_called()
