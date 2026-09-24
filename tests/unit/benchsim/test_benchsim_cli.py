"""scripts/benchsim.py `design` and `generate` subcommands."""

import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "benchsim.py"


def _cli() -> ModuleType:
    spec = spec_from_file_location("benchsim_cli", SCRIPT)
    assert spec and spec.loader
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_requires_subcommand() -> None:
    with pytest.raises(SystemExit) as error:
        _cli().main([])
    assert error.value.code == 2


def test_design_writes_jsonl(tmp_path: Path) -> None:
    rc = _cli().main(
        ["design", "--split", "dev", "--n", "4", "--mutations", "dupC", "--out-root", str(tmp_path)]
    )
    rows = [json.loads(x) for x in (tmp_path / "designs_dev.jsonl").read_text().splitlines()]
    assert rc == 0 and len(rows) == 12 and {r["split"] for r in rows} == {"dev"}


def test_test_split_needs_salt_file_outside_repo(tmp_path: Path) -> None:
    cli = _cli()
    base = ["design", "--split", "test", "--n", "2", "--out-root", str(tmp_path)]
    with pytest.raises(SystemExit):
        cli.main(base)
    inside = SCRIPT.parent / "benchsim.py"  # any file inside the repository
    with pytest.raises(SystemExit):
        cli.main([*base, "--salt-file", str(inside)])
    salt = tmp_path / "salt.txt"
    salt.write_text("secret\n")
    assert cli.main([*base, "--salt-file", str(salt), "--mutations", "dupC"]) == 0
    first = (tmp_path / "designs_test.jsonl").read_text()
    salt.write_text("other\n")
    cli.main([*base, "--salt-file", str(salt), "--mutations", "dupC"])
    assert (tmp_path / "designs_test.jsonl").read_text() != first


def test_generate_writes_manifest_and_fails_on_generation_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli()
    cli.main(
        ["design", "--split", "dev", "--n", "2", "--mutations", "dupC", "--out-root", str(tmp_path)]
    )
    seen: list[Any] = []

    def fake_case(design: Any, ctx: Any) -> dict[str, Any]:
        seen.append(ctx)
        status = "generation_failed" if design.design_id.endswith("0001") else "ok"
        return {"design_id": design.design_id, "status": status}

    monkeypatch.setattr(cli, "require_muconeup", lambda exe: "0.45.0")
    monkeypatch.setattr(cli, "generate_case", fake_case)
    monkeypatch.setattr(cli, "write_variant", lambda *a: (Path("x"), "0" * 64))
    profiles = tmp_path / "profiles_in"
    profiles.mkdir()
    rc = cli.main(
        [
            "generate",
            "--designs",
            str(tmp_path / "designs_dev.jsonl"),
            "--muconeup-config",
            str(tmp_path / "c.json"),
            "--muconeup-profiles",
            str(profiles),
            "--out-root",
            str(tmp_path),
        ]
    )
    manifest = (tmp_path / "dev" / "manifest.jsonl").read_text().splitlines()
    assert rc == 1 and len(manifest) == 6 and len(seen) == 6
    assert seen[0].muconeup_version == "0.45.0" and seen[0].profile_dir == profiles
