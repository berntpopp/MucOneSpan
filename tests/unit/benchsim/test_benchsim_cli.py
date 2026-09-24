"""scripts/benchsim.py `design`, `generate` and `run` subcommands (Git mocked)."""

import json
from dataclasses import asdict
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG

SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "benchsim.py"


def _cli(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """Load the script with a fake worktree ``tmp/wt`` whose main checkout is ``tmp/main``."""
    spec = spec_from_file_location("benchsim_cli", SCRIPT)
    assert spec and spec.loader
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in ("wt", "main", "data"):
        (tmp_path / name).mkdir(exist_ok=True)

    def fake_git(args: list[str], cwd: str | None = None) -> str:
        assert args[:2] == ["git", "rev-parse"]
        if args[2] == "--show-toplevel":
            return f"{tmp_path / 'wt'}\n"
        return f"{tmp_path / 'main' / '.git'}\n"

    monkeypatch.setattr(module, "run_tool", fake_git)
    return module


def test_requires_subcommand(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(SystemExit) as error:
        _cli(tmp_path, monkeypatch).main([])
    assert error.value.code == 2


def test_design_writes_jsonl(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "data"
    rc = _cli(tmp_path, monkeypatch).main(
        ["design", "--split", "dev", "--n", "4", "--mutations", "dupC", "--out-root", str(out)]
    )
    rows = [json.loads(x) for x in (out / "designs_dev_standard.jsonl").read_text().splitlines()]
    assert rc == 0 and len(rows) == 12 and {r["split"] for r in rows} == {"dev"}


def test_design_takes_a_set(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "data"
    cli = _cli(tmp_path, monkeypatch)
    base = ["design", "--split", "dev", "--n", "2", "--mutations", "dupC", "--out-root", str(out)]
    assert cli.main([*base, "--set", "clean"]) == 0
    rows = [json.loads(x) for x in (out / "designs_dev_clean.jsonl").read_text().splitlines()]
    assert {r["bench_set"] for r in rows} == {"clean"}
    assert all(r["design_id"].startswith("dev-clean-") for r in rows)
    with pytest.raises(SystemExit, match="unknown set"):
        cli.main([*base, "--set", "nope"])


def test_bench_config_drives_design(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "data"
    config = tmp_path / "bench.json"
    sizes = {"dev": 2, "val": 2, "test": 2, "stress": 2}
    clean = DEFAULT_BENCH_CONFIG.sets.definitions["clean"]
    profiles = {p: asdict(v) | {"chimera_levels": [0.2]} for p, v in clean.profiles.items()}
    only = {"description": "x", "offpeak_share_cap": None, "profiles": profiles}
    sets = {"definitions": {"only": only}, "default": "only", "headline": "only", "legacy": "only"}
    config.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "design": {"split_sizes": sizes},
                "sets": sets,
                "atlas": {"expected_inconclusive_sets": []},
            }
        )
    )
    args = ["--bench-config", str(config), "design", "--split", "dev", "--mutations", "dupC"]
    rc = _cli(tmp_path, monkeypatch).main([*args, "--out-root", str(out)])
    rows = [json.loads(x) for x in (out / "designs_dev_only.jsonl").read_text().splitlines()]
    assert rc == 0 and len(rows) == 6 and {r["chimera"] for r in rows} == {0.2}
    assert {r["bench_set"] for r in rows} == {"only"}


def test_invalid_bench_config_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = tmp_path / "bench.json"
    config.write_text(json.dumps({"schema_version": 1, "report": {"alpha": 2}}))
    with pytest.raises(SystemExit, match=r"report\.alpha"):
        _cli(tmp_path, monkeypatch).main(["--bench-config", str(config), "preregister"])


def test_bench_config_type_errors_exit_cleanly(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "bench.json"
    config.write_text(json.dumps({"schema_version": 1, "design": {"length_min": [1]}}))
    with pytest.raises(SystemExit, match="length_min"):
        _cli(tmp_path, monkeypatch).main(["--bench-config", str(config), "preregister"])


def test_default_out_root_is_beside_the_repository(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    monkeypatch.chdir(tmp_path / "wt")
    assert cli.main(["design", "--split", "dev", "--n", "1", "--mutations", "dupC"]) == 0
    assert (tmp_path / "MucOneSpan-bench-data" / "designs_dev_standard.jsonl").is_file()


@pytest.mark.parametrize("where", ["wt", "main", "wt/sub"])
@pytest.mark.parametrize("flag", ["--out-root", "--output"])
def test_outputs_inside_repository_or_worktrees_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, where: str, flag: str
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    target = tmp_path / where / ("x.jsonl" if flag == "--output" else "data")
    args = ["design", "--split", "dev", "--n", "1", "--mutations", "dupC", flag, str(target)]
    with pytest.raises(SystemExit, match="outside the repository"):
        cli.main(args)
    assert not target.exists()
    gen = ["generate", "--designs", "d", "--muconeup-config", "c", "--out-root", str(target)]
    if flag == "--out-root":  # refused before designs are read or MucOneUp is called
        with pytest.raises(SystemExit, match="outside the repository"):
            cli.main(gen)


def test_test_split_needs_salt_file_outside_repo(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    out = tmp_path / "data"
    base = ["design", "--split", "test", "--n", "2", "--out-root", str(out)]
    with pytest.raises(SystemExit):
        cli.main(base)
    for where in ("wt", "main"):
        inside = tmp_path / where / "salt.txt"
        inside.write_text("secret\n")
        with pytest.raises(SystemExit, match="outside the repository"):
            cli.main([*base, "--salt-file", str(inside)])
    salt = tmp_path / "salt.txt"
    salt.write_text("secret\n")
    assert cli.main([*base, "--salt-file", str(salt), "--mutations", "dupC"]) == 0
    first = (out / "designs_test_standard.jsonl").read_text()
    salt.write_text("other\n")
    cli.main([*base, "--salt-file", str(salt), "--mutations", "dupC"])
    assert (out / "designs_test_standard.jsonl").read_text() != first
    salt.write_text("\n")
    with pytest.raises(SystemExit, match="empty"):
        cli.main([*base, "--salt-file", str(salt)])


def test_generate_writes_manifest_and_fails_on_generation_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    out = tmp_path / "data"
    cli.main(
        ["design", "--split", "dev", "--n", "2", "--mutations", "dupC", "--out-root", str(out)]
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
            str(out / "designs_dev_standard.jsonl"),
            "--muconeup-config",
            str(tmp_path / "c.json"),
            "--muconeup-profiles",
            str(profiles),
            "--out-root",
            str(out),
        ]
    )
    manifest = (out / "dev" / "manifest.jsonl").read_text().splitlines()
    assert rc == 1 and len(manifest) == 6 and len(seen) == 6
    assert seen[0].muconeup_version == "0.45.0" and seen[0].profile_dir == profiles


def test_generate_stale_case_exits_without_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    out = tmp_path / "data"
    cli.main(
        ["design", "--split", "dev", "--n", "1", "--mutations", "dupC", "--out-root", str(out)]
    )

    def stale(design: Any, ctx: Any) -> dict[str, Any]:
        raise cli.StaleCaseError("case made under other inputs; use a fresh --out-root")

    monkeypatch.setattr(cli, "require_muconeup", lambda exe: "0.45.0")
    monkeypatch.setattr(cli, "generate_case", stale)
    monkeypatch.setattr(cli, "write_variant", lambda *a: (Path("x"), "0" * 64))
    profiles = tmp_path / "profiles_in"
    profiles.mkdir()
    args = ["generate", "--designs", str(out / "designs_dev_standard.jsonl")]
    args += ["--muconeup-config", str(tmp_path / "c.json"), "--muconeup-profiles", str(profiles)]
    with pytest.raises(SystemExit, match="fresh --out-root"):
        cli.main([*args, "--out-root", str(out)])
    assert not (out / "dev" / "manifest.jsonl").exists()


def _manifest(tmp_path: Path) -> Path:
    path = tmp_path / "data" / "dev" / "manifest.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('{"design_id": "dev-a", "profile": "ont_amplicon_r10", "status": "ok"}\n')
    return path


def test_run_derives_engines_results_root_and_models(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    manifest = _manifest(tmp_path)
    calls: list[Any] = []

    def fake_run_split(manifest_arg, engines, results_root, model_for, threads, jobs):
        calls.append((manifest_arg, list(engines), results_root, threads, jobs))
        assert model_for("ont") == "model-ont" and model_for("hifi") == "model-hifi"
        return [{"engine": e, "status": "completed"} for e in engines]

    monkeypatch.setattr(cli, "run_split", fake_run_split)
    rc = cli.main(
        [
            "run",
            "--manifest",
            str(manifest),
            "--model-ont",
            "model-ont",
            "--model-hifi",
            "model-hifi",
        ]
    )
    assert rc == 0 and len(calls) == 1
    _, engines, results_root, threads, jobs = calls[0]
    assert engines == ["ladder"] and threads == 4 and jobs == 1
    assert results_root == tmp_path / "data" / "results" / "dev"


def test_run_accepts_multiple_engines_and_explicit_results_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    manifest = _manifest(tmp_path)
    explicit_root = tmp_path / "data" / "custom_results"
    calls: list[Any] = []
    monkeypatch.setattr(
        cli,
        "run_split",
        lambda manifest_arg, engines, results_root, model_for, threads, jobs: (
            calls.append((list(engines), results_root)) or []
        ),
    )
    rc = cli.main(
        [
            "run",
            "--manifest",
            str(manifest),
            "--engines",
            "ladder,hybrid",
            "--results-root",
            str(explicit_root),
            "--model-ont",
            "m",
            "--model-hifi",
            "m",
        ]
    )
    assert rc == 0
    assert calls == [(["ladder", "hybrid"], explicit_root)]


def test_run_results_root_inside_repository_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    manifest = _manifest(tmp_path)
    with pytest.raises(SystemExit, match="outside the repository"):
        cli.main(
            [
                "run",
                "--manifest",
                str(manifest),
                "--results-root",
                str(tmp_path / "wt" / "res"),
            ]
        )


def test_run_missing_manifest_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli = _cli(tmp_path, monkeypatch)
    with pytest.raises(SystemExit, match="manifest not found"):
        cli.main(["run", "--manifest", str(tmp_path / "missing.jsonl")])


def test_run_without_a_model_fails_only_when_a_case_needs_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    manifest = _manifest(tmp_path)

    def fake_run_split(manifest_arg, engines, results_root, model_for, threads, jobs):
        with pytest.raises(SystemExit, match="model-hifi"):
            model_for("hifi")
        return []

    monkeypatch.setattr(cli, "run_split", fake_run_split)
    assert cli.main(["run", "--manifest", str(manifest), "--model-ont", "m"]) == 0
