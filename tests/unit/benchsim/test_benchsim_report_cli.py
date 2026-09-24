"""scripts/benchsim.py `preregister`, `evaluate`, `realism` and `report` subcommands."""

import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from muc_one_span.benchsim.report import RULE_TEXT
from tests.unit.benchsim.test_benchsim_cli import _cli


def _split(tmp_path: Path, split: str = "dev") -> Path:
    split_dir = tmp_path / "data" / split
    cases = []
    for name, event in (("c1", "dupC"), ("c2", None)):
        case = {
            "design_id": name,
            "profile": "ont_amplicon_r10",
            "status": "ok",
            "design": {
                "design_id": name,
                "profile": "ont_amplicon_r10",
                "event": event,
                "delta_class": "1",
                "depth": 30,
            },
        }
        (split_dir / name).mkdir(parents=True)
        (split_dir / name / "case.json").write_text(json.dumps(case))
        cases.append(case)
    (split_dir / "manifest.jsonl").write_text("".join(json.dumps(c) + "\n" for c in cases))
    return split_dir


def _sample(name: str, truth: str, decision: str, exact: int) -> dict[str, Any]:
    return {
        "sample": name,
        "status": "completed",
        "truth_status": "valid",
        "truth_haplotypes": 2,
        "truth_events": int(truth != "normal"),
        "clinical": {"truth": truth, "decision": decision},
        "alternatives": [
            {
                "pairs": [
                    {"truth": h, "prediction": h, "sequence_exact": bool(exact)}
                    for h in ("h1", "h2")
                ],
                "missing_truth": [],
                "metrics": {"independent_sequence_exact": 2 * exact, "sequence_exact": 2 * exact},
            }
        ],
        "metrics": {
            "all_sequences_exact": {"min": exact, "max": exact},
            "independent_sequence_exact": {"min": 2 * exact, "max": 2 * exact},
            "event_tp": {"min": exact, "max": exact},
            "event_fp": {"min": 0, "max": 0},
        },
    }


def _fake_evaluate(cli: ModuleType, monkeypatch: pytest.MonkeyPatch, calls: list[Any]) -> None:
    def fake_run(args: Any) -> tuple[dict[str, Any], int]:
        calls.append(args)
        engine = args.result_root.name
        exact = 1 if engine == "hybrid" else 0
        return {
            "totals": {},
            "samples": [
                _sample("c1", "pathogenic", "PATHOGENIC", exact),
                _sample("c2", "normal", "NO_PATHOGENIC_VARIANT_DETECTED", exact),
            ],
        }, 0

    monkeypatch.setattr(cli, "evaluate_run", fake_run)


def test_preregister_writes_ledger_under_out_root_test(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    assert cli.main(["preregister", "--out-root", str(tmp_path / "data")]) == 0
    entry = json.loads((tmp_path / "data" / "test" / "preregistration.jsonl").read_text())
    assert entry["rule_text"] == RULE_TEXT


def test_preregister_inside_repository_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    with pytest.raises(SystemExit, match="outside the repository"):
        _cli(tmp_path, monkeypatch).main(["preregister", "--out-root", str(tmp_path / "wt")])


def test_evaluate_test_split_refused_before_reading_truth(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    calls: list[Any] = []
    _fake_evaluate(cli, monkeypatch, calls)
    with pytest.raises(SystemExit, match="pre-regist"):
        cli.main(["evaluate", "--split", "test", "--out-root", str(tmp_path / "data")])
    assert calls == []
    _split(tmp_path, "test")
    for engine in ("ladder",):
        (tmp_path / "data" / "results" / "test" / engine).mkdir(parents=True)
    cli.main(["preregister", "--out-root", str(tmp_path / "data")])
    assert cli.main(["evaluate", "--split", "test", "--out-root", str(tmp_path / "data")]) == 0
    assert len(calls) == 1
    evaluation = json.loads(
        (tmp_path / "data" / "results" / "test" / "ladder" / "evaluation.json").read_text()
    )
    audit = evaluation["preregistration"]
    assert set(audit) == {"sha256", "registered_at", "test_first_evaluated_at"}
    first = audit["test_first_evaluated_at"]
    # A second evaluation keeps the first unsealing time; re-registration is refused.
    cli.main(["evaluate", "--split", "test", "--out-root", str(tmp_path / "data")])
    again = json.loads(
        (tmp_path / "data" / "results" / "test" / "ladder" / "evaluation.json").read_text()
    )
    assert again["preregistration"]["test_first_evaluated_at"] == first
    with pytest.raises(SystemExit, match="already evaluated"):
        cli.main(["preregister", "--out-root", str(tmp_path / "data")])
    assert cli.main(["report", "--split", "test", "--out-root", str(tmp_path / "data")]) == 0
    report = json.loads((tmp_path / "data" / "results" / "test" / "report.json").read_text())
    assert report["preregistration"] == audit


def test_evaluate_loads_scripts_evaluate_lazily(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    assert cli.evaluate_run is None and callable(cli._load_evaluate().run)
    _split(tmp_path)
    (tmp_path / "data" / "results" / "dev" / "ladder").mkdir(parents=True)
    loaded: list[Any] = []

    class FakeModule:
        @staticmethod
        def run(args: Any) -> tuple[dict[str, Any], int]:
            loaded.append(args)
            return {"samples": []}, 1

    monkeypatch.setattr(cli, "_load_evaluate", lambda: FakeModule)
    assert cli.main(["evaluate", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 1
    assert len(loaded) == 1


def test_evaluate_then_report_writes_json_and_markdown(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    calls: list[Any] = []
    _fake_evaluate(cli, monkeypatch, calls)
    split_dir = _split(tmp_path)
    results = tmp_path / "data" / "results" / "dev"
    for engine in ("ladder", "hybrid"):
        (results / engine).mkdir(parents=True)
    rc = cli.main(
        [
            "evaluate",
            "--split",
            "dev",
            "--engines",
            "ladder,hybrid",
            "--out-root",
            str(tmp_path / "data"),
        ]
    )
    assert rc == 0 and [c.result_root.name for c in calls] == ["ladder", "hybrid"]
    assert calls[0].expected_samples == results / "ladder" / "inventory_dev.json"
    assert calls[0].truth_root == split_dir
    assert (results / "hybrid" / "evaluation.json").is_file()
    rc = cli.main(
        [
            "report",
            "--split",
            "dev",
            "--baseline",
            "ladder",
            "--candidate",
            "hybrid",
            "--out-root",
            str(tmp_path / "data"),
        ]
    )
    report = json.loads((results / "report.json").read_text())
    assert rc == 0 and report["decision"]["adopt"] is False  # 2 cases cannot reach significance
    assert report["engines"]["hybrid"]["rows"][0]["alleles"][0]["allele_exact"] == 1
    assert report["decision"]["profiles"]["ont_amplicon_r10"]["n_alleles"] == 4
    tables = report["tables"]["hybrid"]
    assert set(tables) == {"stratified", "pooled", "events", "confusion"}
    assert tables["pooled"]["ont_amplicon_r10"]["allele_exact"]["point"] == 1.0
    assert report["preregistration"] is None
    text = (results / "report.md").read_text()
    assert "NOT ADOPTED" in text and "failure atlas: composition" in text


def test_report_invalid_truth_is_a_visible_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _split(tmp_path)
    engine = tmp_path / "data" / "results" / "dev" / "ladder"
    engine.mkdir(parents=True)
    bad = {"sample": "c1", "status": "invalid_truth", "truth_status": "invalid"}
    (engine / "evaluation.json").write_text(json.dumps({"samples": [bad]}))
    with pytest.raises(SystemExit, match="no valid truth"):
        _cli(tmp_path, monkeypatch).main(
            ["report", "--split", "dev", "--out-root", str(tmp_path / "data")]
        )


def test_report_single_engine_has_no_decision(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _fake_evaluate(cli, monkeypatch, [])
    _split(tmp_path)
    (tmp_path / "data" / "results" / "dev" / "ladder").mkdir(parents=True)
    cli.main(["evaluate", "--split", "dev", "--out-root", str(tmp_path / "data")])
    assert cli.main(["report", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 0
    report = json.loads((tmp_path / "data" / "results" / "dev" / "report.json").read_text())
    assert report["decision"] is None and "ladder" in report["engines"]


def test_report_missing_evaluation_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _split(tmp_path)
    with pytest.raises(SystemExit, match=r"evaluation\.json"):
        _cli(tmp_path, monkeypatch).main(
            ["report", "--split", "dev", "--out-root", str(tmp_path / "data")]
        )


def test_realism_aggregates_per_profile_and_keeps_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _split(tmp_path)
    seen: list[Path] = []

    def fake_case_metrics(
        case_dir: Path, muconeup_config: Any = None, flank_fasta: Any = None
    ) -> dict[str, Any]:
        seen.append(case_dir)
        if case_dir.name == "c2":
            raise ValueError("broken truth")
        return {"n_reads": 3}

    monkeypatch.setattr(cli, "case_metrics", fake_case_metrics)
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases: {"n_cases": len(cases)})
    monkeypatch.setattr(cli, "realism_compare", lambda m, t, p: {"check": {"pass": True}})
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    assert cli.main(["realism", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 0
    out = json.loads((tmp_path / "data" / "dev" / "realism.json").read_text())
    prof = out["profiles"]["ont_amplicon_r10"]
    assert prof["aggregate"] == {"n_cases": 1} and prof["compare"]["check"]["pass"] is True
    assert out["failures"] == [{"design_id": "c2", "error": "ValueError: broken truth"}]
    assert "ont_amplicon_r10" in (tmp_path / "data" / "dev" / "realism.md").read_text()


def test_realism_profile_without_targets(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _split(tmp_path)
    monkeypatch.setattr(cli, "case_metrics", lambda *a, **k: {})
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases: {"n_cases": len(cases)})

    def no_section(m: Any, t: Any, p: str) -> dict[str, Any]:
        raise KeyError(p)

    monkeypatch.setattr(cli, "realism_compare", no_section)
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    assert cli.main(["realism", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 0
    out = json.loads((tmp_path / "data" / "dev" / "realism.json").read_text())
    assert out["profiles"]["ont_amplicon_r10"]["compare"] is None
