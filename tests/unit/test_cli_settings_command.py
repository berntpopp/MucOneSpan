"""Tests for the ``muconespan settings show|validate`` commands."""

from __future__ import annotations

import json
from dataclasses import fields, is_dataclass, replace
from pathlib import Path

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings, load_settings

EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "runtime-settings.json"
SECTIONS = [
    f.name for f in fields(RuntimeSettings) if is_dataclass(getattr(DEFAULT_SETTINGS, f.name))
]


def _invoke(*args: str) -> tuple[int, str, str]:
    result = CliRunner().invoke(main, list(args))
    return result.exit_code, result.stdout, result.stderr


def test_show_without_config_prints_defaults_that_load_back(tmp_path: Path) -> None:
    code, out, _ = _invoke("settings", "show")
    assert code == 0
    path = tmp_path / "shown.json"
    path.write_text(out)
    assert json.loads(out)["schema_version"] == DEFAULT_SETTINGS.schema_version
    assert load_settings(path) == DEFAULT_SETTINGS


def test_example_configuration_is_in_sync_with_settings_show() -> None:
    code, out, _ = _invoke("settings", "show")
    assert code == 0
    assert EXAMPLE.read_text() == out


def test_show_with_config_prints_effective_settings(tmp_path: Path) -> None:
    threads = DEFAULT_SETTINGS.run.threads + 1
    seed = DEFAULT_SETTINGS.hybrid.seed + 1
    config = tmp_path / "in.json"
    config.write_text(
        json.dumps({"schema_version": 1, "run": {"threads": threads}, "hybrid": {"seed": seed}})
    )
    code, out, _ = _invoke("settings", "show", "--config", str(config))
    assert code == 0
    shown = tmp_path / "shown.json"
    shown.write_text(out)
    expected = replace(
        DEFAULT_SETTINGS,
        run=replace(DEFAULT_SETTINGS.run, threads=threads),
        hybrid=replace(DEFAULT_SETTINGS.hybrid, seed=seed),
    )
    assert load_settings(shown) == expected


def test_show_uses_global_config_when_no_local_config(tmp_path: Path) -> None:
    threads = DEFAULT_SETTINGS.run.threads + 1
    config = tmp_path / "in.json"
    config.write_text(json.dumps({"schema_version": 1, "run": {"threads": threads}}))
    code, out, _ = _invoke("--config", str(config), "settings", "show")
    assert code == 0
    assert json.loads(out)["run"]["threads"] == threads


@pytest.mark.parametrize("section", SECTIONS)
def test_show_section_prints_one_loadable_section(tmp_path: Path, section: str) -> None:
    code, out, _ = _invoke("settings", "show", "--section", section)
    assert code == 0
    data = json.loads(out)
    assert set(data) == {"schema_version", section}
    path = tmp_path / "section.json"
    path.write_text(out)
    assert getattr(load_settings(path), section) == getattr(DEFAULT_SETTINGS, section)


def test_show_rejects_unknown_section() -> None:
    code, _, _ = _invoke("settings", "show", "--section", "not_a_section")
    assert code != 0


def test_show_reports_invalid_config(tmp_path: Path) -> None:
    config = tmp_path / "bad.json"
    config.write_text(json.dumps({"schema_version": 1, "run": {"threads": 0}}))
    code, out, err = _invoke("settings", "show", "--config", str(config))
    assert code != 0
    assert out == ""
    assert "threads" in err


def test_validate_accepts_valid_file() -> None:
    code, out, _ = _invoke("settings", "validate", str(EXAMPLE))
    assert code == 0
    assert "valid" in out


@pytest.mark.parametrize(
    ("contents", "fragment"),
    [
        ('{"schema_version": 1, "run": {"bogus": 1}}', "Unknown run fields: bogus"),
        ('{"run": {}}', "schema_version"),
        ('{"schema_version": 1, "hybrid": {"seed": "x"}}', "seed"),
        ('{"schema_version": 1, "schema_version": 1}', "Duplicate configuration key"),
        ("[1]", "JSON object"),
        ("{not json", "Expecting property name"),
    ],
)
def test_validate_reports_first_error_and_fails(
    tmp_path: Path, contents: str, fragment: str
) -> None:
    config = tmp_path / "bad.json"
    config.write_text(contents)
    code, out, err = _invoke("settings", "validate", str(config))
    assert code != 0
    assert fragment in err
    assert str(config) in err
    assert out == ""


def test_validate_reports_missing_file(tmp_path: Path) -> None:
    code, _, err = _invoke("settings", "validate", str(tmp_path / "missing.json"))
    assert code != 0
    assert "missing.json" in err
