"""scripts/benchsim.py ``calibrate --stage lengths`` (Task 15d).

Unlike `test_benchsim_calibrate_cli.py`, the engine run itself is not mocked here:
these tests build tiny real cases (real bundled repeat dictionary, synthetic
spanning reads) and run the real hybrid length model end to end through
`benchsim calibrate --stage lengths` and `benchsim calibrate-report`, proving the
new stage plugs into the unchanged 15b/15c grid/point/manifest/resume machinery
and the 15c report path exactly like a full-pipeline calibration.
"""

import json
from pathlib import Path

import pytest

from muc_one_span.settings import load_settings
from tests.unit.benchsim.test_benchsim_cli import _cli
from tests.unit.benchsim.test_calibration_lengths import _case, _manifest

BASE_ARGS = ("--engine", "hybrid", "--stage", "lengths")


def _grid(tmp_path: Path, data: dict) -> Path:
    path = tmp_path / "grid.json"
    path.write_text(json.dumps(data))
    return path


def _one_case_split(tmp_path: Path) -> Path:
    split_dir = tmp_path / "data" / "dev"
    rows = [_case(split_dir, "c1", [["X"] * 30], 150, seed=1)]
    _manifest(split_dir, rows)
    return split_dir


def test_stage_lengths_requires_engine_hybrid(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    grid = _grid(tmp_path, {})
    with pytest.raises(SystemExit, match="requires --engine hybrid"):
        cli.main(
            [
                "calibrate",
                "--split",
                "dev",
                "--engine",
                "ladder",
                "--stage",
                "lengths",
                "--grid",
                str(grid),
                "--out-root",
                str(tmp_path / "data"),
            ]
        )
    assert not (tmp_path / "data" / "calibration").exists()


def test_stage_lengths_refuses_non_length_grid_keys_before_any_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    (tmp_path / "data" / "dev").mkdir(parents=True)
    (tmp_path / "data" / "dev" / "manifest.jsonl").write_text("")
    grid = _grid(tmp_path, {"hybrid.n_poa": [10]})  # not a spans/lengths/smear setting
    with pytest.raises(SystemExit, match=r"not a length-model key.*hybrid\.n_poa"):
        cli.main(
            [
                "calibrate",
                "--split",
                "dev",
                *BASE_ARGS,
                "--grid",
                str(grid),
                "--out-root",
                str(tmp_path / "data"),
            ]
        )
    assert not (tmp_path / "data" / "calibration").exists()


def test_calibrate_stage_lengths_runs_the_real_length_model_and_reports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _one_case_split(tmp_path)
    grid = _grid(tmp_path, {"hybrid.min_peak_reads": [8, 4]})
    base = ["--out-root", str(tmp_path / "data"), *BASE_ARGS, "--grid", str(grid)]
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0
    cal_dir = tmp_path / "data" / "calibration" / "dev" / "grid"
    manifest = json.loads((cal_dir / "calibration.json").read_text())
    assert manifest["inputs"]["stage"] == "lengths"
    assert [p["status"] for p in manifest["points"]] == ["evaluated", "evaluated"]
    for point in manifest["points"]:
        evaluation = json.loads(
            (cal_dir / point["sha256"] / "results" / "hybrid" / "evaluation.json").read_text()
        )
        (row,) = evaluation["rows"]
        assert row["sample"] == "c1" and row["status"] == "ok"
        assert row["truth_allele_count"] == 1 and row["allele_count_exact"] == 1

    objective = tmp_path / "objective.json"
    objective.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "constraints": {"missed_alleles": {"max": 0}},
                "rank": ["-allele_length_exact", "case_length_exact"],
            }
        )
    )
    report_args = [
        "calibrate-report",
        "--split",
        "dev",
        "--name",
        "grid",
        "--objective",
        str(objective),
        "--out-root",
        str(tmp_path / "data"),
    ]
    assert cli.main(report_args) == 0
    report = json.loads((cal_dir / "calibration-report.json").read_text())
    assert report["stage"] == "lengths"
    assert report["recommended"] is not None
    assert "allele_length_exact" in (cal_dir / "calibration-report.md").read_text()
    recommended = load_settings(cal_dir / "recommended-config.json")
    assert recommended.hybrid.min_peak_reads in (4, 8)
    assert recommended.run.engine == "hybrid"


def test_calibrate_stage_lengths_resumes_without_rerunning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _one_case_split(tmp_path)
    calls: list[int] = []
    real = cli.run_lengths_stage

    def counted(*args: object, **kwargs: object) -> object:
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "run_lengths_stage", counted)
    grid = _grid(tmp_path, {})
    base = ["--out-root", str(tmp_path / "data"), *BASE_ARGS, "--grid", str(grid)]
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0
    assert len(calls) == 1
    assert cli.main(["calibrate", "--split", "dev", *base]) == 0
    assert len(calls) == 1  # the one point was already evaluated: not rerun


def test_calibrate_refuses_a_stage_change_under_the_same_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _one_case_split(tmp_path)
    grid = _grid(tmp_path, {})
    base = ["--out-root", str(tmp_path / "data"), "--engine", "hybrid", "--grid", str(grid)]
    assert cli.main(["calibrate", "--split", "dev", "--stage", "lengths", *base]) == 0
    with pytest.raises(SystemExit, match="stage"):
        cli.main(["calibrate", "--split", "dev", *base])  # default --stage full
