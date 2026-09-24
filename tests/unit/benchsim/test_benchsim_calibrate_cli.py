"""scripts/benchsim.py `calibrate` and `calibrate-report` (engine runs and scoring mocked)."""

import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from muc_one_span.settings import DEFAULT_SETTINGS, load_settings
from tests.unit.benchsim.test_benchsim_cli import _cli
from tests.unit.benchsim.test_benchsim_report_cli import _sample, _split

HET = DEFAULT_SETTINGS.hybrid.het_af_min
BETTER = HET / 2  # the synthetic evaluator scores this value exact, the default not
GRID = {"hybrid.het_af_min": [HET, BETTER]}
OBJECTIVE = {
    "schema_version": 1,
    "constraints": {"clinical_false_negative": {"max": 0}},
    "rank": ["-per_allele_exact", "inconclusive_rate"],
}


def _json(path: Path, data: object) -> Path:
    path.write_text(json.dumps(data))
    return path


def _mock_engine(cli: ModuleType, monkeypatch: pytest.MonkeyPatch) -> dict[str, list[Any]]:
    calls: dict[str, list[Any]] = {"run": [], "evaluate": []}

    def fake_run_split(manifest, engines, results_root, model_for, threads, jobs, config=None):
        calls["run"].append((engines, results_root, config))
        (results_root / engines[0]).mkdir(parents=True, exist_ok=True)
        return []

    def fake_evaluate(args: Any) -> tuple[dict[str, Any], int]:
        calls["evaluate"].append(args)
        config = json.loads((args.result_root.parents[1] / "config.json").read_text())
        exact = int(config["hybrid"]["het_af_min"] == BETTER)
        cfn_decision = "PATHOGENIC" if exact else "NO_PATHOGENIC_VARIANT_DETECTED"
        samples = [
            _sample("c1", "pathogenic", cfn_decision, exact),
            _sample("c2", "normal", "NO_PATHOGENIC_VARIANT_DETECTED", exact),
        ]
        return {"totals": {}, "samples": samples}, 0

    monkeypatch.setattr(cli, "run_split", fake_run_split)
    monkeypatch.setattr(cli, "evaluate_run", fake_evaluate)
    return calls


def _setup(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[ModuleType, dict, list]:
    cli = _cli(tmp_path, monkeypatch)
    for split in ("dev", "val"):
        _split(tmp_path, split)
    calls = _mock_engine(cli, monkeypatch)
    grid = _json(tmp_path / "grid.json", GRID)
    base = ["--out-root", str(tmp_path / "data"), "--engine", "hybrid", "--grid", str(grid)]
    return cli, calls, base


def _cal_dir(tmp_path: Path, split: str = "dev", name: str = "grid") -> Path:
    return tmp_path / "data" / "calibration" / split / name


def test_calibrate_runs_every_point_through_run_and_evaluate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, calls, base = _setup(tmp_path, monkeypatch)
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0
    manifest = json.loads((_cal_dir(tmp_path) / "calibration.json").read_text())
    inputs = manifest["inputs"]
    assert inputs["split"] == "dev" and inputs["engine"] == "hybrid"
    assert inputs["grid"] == GRID and inputs["base_config"] is None
    assert set(inputs["versions"]) >= {"muconespan", "muconeup"}
    assert set(inputs["design_seeds"]) == {"c1", "c2"}
    assert [p["status"] for p in manifest["points"]] == ["evaluated", "evaluated"]
    assert len(calls["run"]) == len(calls["evaluate"]) == len(GRID["hybrid.het_af_min"])
    for point in manifest["points"]:
        point_dir = _cal_dir(tmp_path) / point["sha256"]
        loaded = load_settings(point_dir / "config.json")
        assert loaded.hybrid.het_af_min == point["values"]["hybrid.het_af_min"]
        assert loaded.run.engine == "hybrid"
        assert (point_dir / "results" / "hybrid" / "evaluation.json").is_file()
    engines, results_root, config = calls["run"][0]
    assert engines == ["hybrid"] and config == results_root.parent / "config.json"


def test_calibrate_resumes_without_rerunning_completed_points(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, calls, base = _setup(tmp_path, monkeypatch)
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0
    first = len(calls["run"])
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0
    assert len(calls["run"]) == first


def test_calibrate_reruns_a_failed_point(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli, calls, base = _setup(tmp_path, monkeypatch)
    real = cli.run_split

    def flaky(*args: Any, **kwargs: Any) -> list[Any]:
        if not calls["run"]:
            calls["run"].append("boom")
            raise RuntimeError("disk full")
        return list(real(*args, **kwargs))

    monkeypatch.setattr(cli, "run_split", flaky)
    assert cli.main(["calibrate", "--split", "dev", *base]) == 1
    manifest = json.loads((_cal_dir(tmp_path) / "calibration.json").read_text())
    assert [p["status"] for p in manifest["points"]] == ["failed", "evaluated"]
    assert "disk full" in manifest["points"][0]["error"]
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0


def test_calibrate_refuses_changed_inputs_under_the_same_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0
    with pytest.raises(SystemExit, match="engine"):
        cli.main(["calibrate", "--split", "dev", *base[:2], "--engine", "ladder", *base[4:]])


def test_calibrate_refuses_the_sealed_test_split(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, calls, base = _setup(tmp_path, monkeypatch)
    _split(tmp_path, "test")
    with pytest.raises(SystemExit, match="sealed"):
        cli.main(["calibrate", "--split", "test", *base])
    assert calls["run"] == []
    with pytest.raises(SystemExit, match="sealed"):
        cli.main(["calibrate-report", "--split", "test", "--name", "grid", "--objective", "x"])
    with pytest.raises(SystemExit, match="not 'stress'"):
        cli.main(["calibrate", "--split", "stress", *base])


def test_invalid_grid_fails_before_any_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli, calls, base = _setup(tmp_path, monkeypatch)
    _json(tmp_path / "grid.json", {"hybrid.het_af_min": [HET, -1.0]})
    with pytest.raises(SystemExit, match=r"hybrid\.het_af_min"):
        cli.main(["calibrate", "--split", "dev", *base])
    _json(tmp_path / "grid.json", {"hybrid.nope": [1]})
    with pytest.raises(SystemExit, match="Unknown hybrid fields"):
        cli.main(["calibrate", "--split", "dev", *base])
    assert calls["run"] == []
    assert not _cal_dir(tmp_path).exists()


def test_calibrate_rejects_unsafe_names_and_a_missing_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    with pytest.raises(SystemExit, match="plain directory name"):
        cli.main(["calibrate", "--split", "dev", *base, "--name", "../x"])
    (tmp_path / "data" / "val" / "manifest.jsonl").unlink()
    with pytest.raises(SystemExit, match="manifest not found"):
        cli.main(["calibrate", "--split", "val", *base])


def _report(cli: ModuleType, tmp_path: Path, *extra: str, split: str = "dev") -> int:
    objective = _json(tmp_path / "objective.json", OBJECTIVE)
    args = ["calibrate-report", "--split", split, "--out-root", str(tmp_path / "data")]
    code: int = cli.main([*args, "--name", "grid", "--objective", str(objective), *extra])
    return code


def test_calibrate_report_ranks_and_recommends_a_loadable_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    cli.main(["calibrate", "--split", "dev", *base])
    assert _report(cli, tmp_path) == 0
    cal = _cal_dir(tmp_path)
    report = json.loads((cal / "calibration-report.json").read_text())
    best, worst = report["points"]
    assert best["values"] == {"hybrid.het_af_min": BETTER} and best["feasible"]
    assert worst["violations"] == ["clinical_false_negative: value 1 > max 0"]
    assert report["recommended"] == best["sha256"]
    metric = best["metrics"]["per_allele_exact"]
    assert metric["ci_low"] <= metric["value"] <= metric["ci_high"]
    recommended = load_settings(cal / "recommended-config.json")
    assert recommended.hybrid.het_af_min == BETTER
    assert (cal / "recommended-config.json").read_text() == (
        cal / best["sha256"] / "config.json"
    ).read_text()
    provenance = json.loads((cal / "recommended-config.provenance.json").read_text())
    assert provenance["point_sha256"] == best["sha256"]
    assert provenance["objective"] == OBJECTIVE and provenance["grid"] == GRID
    assert provenance["split"] == "dev"
    grid = json.loads((cal / "recommended-grid.json").read_text())
    assert grid == {"hybrid.het_af_min": [BETTER]}
    markdown = (cal / "calibration-report.md").read_text()
    assert best["sha256"] in markdown and "per_allele_exact" in markdown


def test_confirmation_on_val_shows_the_dev_to_val_shift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    cli.main(["calibrate", "--split", "dev", *base])
    _report(cli, tmp_path)
    grid = str(_cal_dir(tmp_path) / "recommended-grid.json")
    confirm = [*base[:4], "--grid", grid, "--name", "grid"]
    assert cli.main(["calibrate", "--split", "val", *confirm]) == 0
    assert _report(cli, tmp_path, "--shift-from", "grid", split="val") == 0
    report = json.loads((_cal_dir(tmp_path, "val") / "calibration-report.json").read_text())
    (point,) = report["points"]
    shift = report["shift"]["points"][point["sha256"]]
    assert shift["per_allele_exact"]["dev"] == shift["per_allele_exact"]["val"]
    assert shift["per_allele_exact"]["delta"] == 0
    assert "dev" in (_cal_dir(tmp_path, "val") / "calibration-report.md").read_text()


def test_shift_is_only_for_the_confirmation_split(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    cli.main(["calibrate", "--split", "dev", *base])
    with pytest.raises(SystemExit, match="--shift-from"):
        _report(cli, tmp_path, "--shift-from", "grid")


def test_no_feasible_point_writes_no_recommendation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    _json(tmp_path / "grid.json", {"hybrid.het_af_min": [HET]})
    cli.main(["calibrate", "--split", "dev", *base])
    assert _report(cli, tmp_path) == 1
    cal = _cal_dir(tmp_path)
    assert json.loads((cal / "calibration-report.json").read_text())["recommended"] is None
    assert not (cal / "recommended-config.json").exists()


def test_calibrate_report_needs_a_calibration_and_a_valid_objective(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    with pytest.raises(SystemExit, match=r"calibration\.json"):
        _report(cli, tmp_path)
    cli.main(["calibrate", "--split", "dev", *base])
    objective = _json(tmp_path / "bad-objective.json", {"schema_version": 1})
    args = ["calibrate-report", "--split", "dev", "--out-root", str(tmp_path / "data")]
    with pytest.raises(SystemExit, match="rank"):
        cli.main([*args, "--name", "grid", "--objective", str(objective)])


def test_report_lists_failed_points_and_unmatched_shift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, calls, base = _setup(tmp_path, monkeypatch)
    cli.main(["calibrate", "--split", "dev", *base])
    real = cli.run_split

    def fail_other(manifest, engines, results_root, *args: Any, **kwargs: Any) -> list[Any]:
        config = json.loads(kwargs["config"].read_text())
        if config["hybrid"]["het_af_min"] != BETTER:
            raise OSError("no space")
        return list(real(manifest, engines, results_root, *args, **kwargs))

    monkeypatch.setattr(cli, "run_split", fail_other)
    _json(tmp_path / "grid.json", {"hybrid.het_af_min": [BETTER, HET / 4]})
    assert cli.main(["calibrate", "--split", "val", *base]) == 1
    assert _report(cli, tmp_path, "--shift-from", "grid", split="val") == 0
    cal = _cal_dir(tmp_path, "val")
    report = json.loads((cal / "calibration-report.json").read_text())
    (failed,) = report["not_evaluated"]
    assert failed["status"] == "failed" and "no space" in failed["error"]
    markdown = (cal / "calibration-report.md").read_text()
    assert "## Not evaluated" in markdown
    assert calls["run"]


def test_report_shift_without_a_matching_dev_point(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    cli.main(["calibrate", "--split", "dev", *base])
    _json(tmp_path / "grid.json", {"hybrid.het_af_min": [BETTER / 2]})
    cli.main(["calibrate", "--split", "val", *base])
    _report(cli, tmp_path, "--shift-from", "grid", split="val")
    cal = _cal_dir(tmp_path, "val")
    shift = json.loads((cal / "calibration-report.json").read_text())["shift"]
    assert list(shift["points"].values()) == [None]
    assert "no dev point" in (cal / "calibration-report.md").read_text()


def test_report_refuses_rows_it_cannot_normalize(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    cli.main(["calibrate", "--split", "dev", *base])
    (tmp_path / "data" / "dev" / "manifest.jsonl").write_text("")
    with pytest.raises(SystemExit, match="cannot normalize"):
        _report(cli, tmp_path)


def test_a_tampered_point_config_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, _, base = _setup(tmp_path, monkeypatch)
    cli.main(["calibrate", "--split", "dev", *base])
    manifest = json.loads((_cal_dir(tmp_path) / "calibration.json").read_text())
    point = _cal_dir(tmp_path) / manifest["points"][0]["sha256"]
    (point / "results" / "hybrid" / "evaluation.json").unlink()
    (point / "config.json").write_text("{}")
    with pytest.raises(SystemExit, match="content address"):
        cli.main(["calibrate", "--split", "dev", *base])
