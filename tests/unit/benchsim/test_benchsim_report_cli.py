"""scripts/benchsim.py `preregister`, `evaluate`, `realism` and `report` subcommands."""

import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG
from muc_one_span.benchsim.report import RULE_TEXT, rule_sha256
from tests.unit.benchsim.test_benchsim_cli import _cli

STANDARD = DEFAULT_BENCH_CONFIG.sets.headline
FIXTURE = (("c1", "dupC", STANDARD), ("c2", None, STANDARD))


def _split(
    tmp_path: Path, split: str = "dev", fixture: tuple[tuple[Any, ...], ...] = FIXTURE
) -> Path:
    split_dir = tmp_path / "data" / split
    cases = []
    for name, event, bench_set in fixture:
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
                **({"bench_set": bench_set} if bench_set else {}),
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
    assert report["headline_set"] == STANDARD and report["set_order"] == [STANDARD]
    tables = report["sets"][STANDARD]["tables"]["hybrid"]
    assert set(tables) == {"stratified", "pooled", "events", "confusion"}
    assert tables["pooled"]["ont_amplicon_r10"]["allele_exact"]["point"] == 1.0
    assert report["preregistration"] is None
    # Task 12e: the absolute targets are evaluated per set, not restricted to `standard`
    # (this fixture has no `clean` cases at all: not present, never FAIL; blocks adoption).
    targets = report["decision"]["targets"]
    assert set(targets) == {STANDARD, "clean"}
    assert targets[STANDARD]["pass"] is True
    assert targets["clean"]["present"] is False and targets["clean"]["pass"] is None
    text = (results / "report.md").read_text()
    assert "NOT ADOPTED" in text and "failure atlas: composition" in text
    assert "Part 2: absolute targets" in text and "Set verdict:" in text


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
    text = (tmp_path / "data" / "results" / "dev" / "report.md").read_text()
    assert "ADOPT" not in text and "no candidate" in text


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
        case_dir: Path, muconeup_config: Any = None, flank_fasta: Any = None, cfg: Any = None
    ) -> dict[str, Any]:
        seen.append(case_dir)
        if case_dir.name == "c2":
            raise ValueError("broken truth")
        return {"n_reads": 3}

    monkeypatch.setattr(cli, "case_metrics", fake_case_metrics)
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases, cfg: {"n_cases": len(cases)})
    monkeypatch.setattr(cli, "realism_compare", lambda m, t, p, cfg: {"check": {"pass": True}})
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    assert cli.main(["realism", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 0
    out = json.loads((tmp_path / "data" / "dev" / "realism.json").read_text())
    prof = out["sets"][STANDARD]["profiles"]["ont_amplicon_r10"]
    assert prof["aggregate"] == {"n_cases": 1} and prof["compare"]["check"]["pass"] is True
    assert out["failures"] == [{"design_id": "c2", "error": "ValueError: broken truth"}]
    assert out["bench_config_sha256"] == DEFAULT_BENCH_CONFIG.sha256()
    text = (tmp_path / "data" / "dev" / "realism.md").read_text()
    assert "ont_amplicon_r10" in text and "indicative" in text.lower()


def test_realism_profile_without_targets(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _split(tmp_path)
    monkeypatch.setattr(cli, "case_metrics", lambda *a, **k: {})
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases, cfg: {"n_cases": len(cases)})

    def no_section(m: Any, t: Any, p: str, cfg: Any) -> dict[str, Any]:
        raise KeyError(p)

    monkeypatch.setattr(cli, "realism_compare", no_section)
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    assert cli.main(["realism", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 0
    out = json.loads((tmp_path / "data" / "dev" / "realism.json").read_text())
    assert out["sets"][STANDARD]["profiles"]["ont_amplicon_r10"]["compare"] is None


def test_realism_test_split_needs_preregistration_and_marks_first_evaluation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _split(tmp_path, "test")
    monkeypatch.setattr(cli, "case_metrics", lambda *a, **k: {})
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases, cfg: {"n_cases": len(cases)})
    monkeypatch.setattr(cli, "realism_compare", lambda m, t, p, cfg: {})
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    argv = ["realism", "--split", "test", "--out-root", str(tmp_path / "data")]
    with pytest.raises(SystemExit, match="pre-regist"):
        cli.main(argv)
    marker = tmp_path / "data" / "test" / "first_evaluation.json"
    assert not marker.exists()
    cli.main(["preregister", "--out-root", str(tmp_path / "data")])
    assert cli.main(argv) == 0
    assert marker.is_file()  # realism reads test truth, so it unseals test too
    with pytest.raises(SystemExit, match="already evaluated"):
        cli.main(["preregister", "--out-root", str(tmp_path / "data")])


def test_realism_settings_mismatch_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _split(tmp_path)
    monkeypatch.setattr(cli, "case_metrics", lambda *a, **k: {})
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases, cfg: {"n_cases": len(cases)})

    def mismatch(m: Any, t: Any, p: str, cfg: Any) -> dict[str, Any]:
        raise ValueError("histogram keys ['lt60u', 'ge60u']")

    monkeypatch.setattr(cli, "realism_compare", mismatch)
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    with pytest.raises(SystemExit, match="lt60u"):
        cli.main(["realism", "--split", "dev", "--out-root", str(tmp_path / "data")])


def test_report_writes_the_reason_atlas_with_config_thresholds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    reason = "Allele 1: 7 primary alignments, below the per-allele depth gate (30)."

    def fake_run(args: Any) -> tuple[dict[str, Any], int]:
        unclear = _sample("c1", "pathogenic", "INCONCLUSIVE", 0)
        unclear["clinical"]["reasons"] = [reason]
        unclear["reconstruction_flags"] = ["ambiguous_reconstruction", "iupac_bases"]
        return {"samples": [unclear, _sample("c2", "normal", "INCONCLUSIVE", 0)]}, 0

    monkeypatch.setattr(cli, "evaluate_run", fake_run)
    _split(tmp_path)
    (tmp_path / "data" / "results" / "dev" / "ladder").mkdir(parents=True)
    data = str(tmp_path / "data")
    cli.main(["evaluate", "--split", "dev", "--out-root", data])
    assert cli.main(["report", "--split", "dev", "--out-root", data]) == 0
    report = json.loads((tmp_path / "data" / "results" / "dev" / "report.json").read_text())
    atlas = report["sets"][STANDARD]["atlas"]["ladder"]
    assert (atlas["n_cases"], atlas["n_atlas"]) == (2, 2)
    keys = {r["reason"]: r["k"] for r in atlas["reasons"]}
    assert keys["evaluator: iupac_bases"] == 1 and keys["unrecorded"] == 1  # legacy c2
    assert keys["gate: # primary alignments, below the per-allele depth gate (#)"] == 1
    total = atlas["split_summary"][-1]  # the fixture records no realized depth
    assert (total["expected"], total["resolvable"], total["depth_unknown"]) == (0, 0, 2)
    text = (tmp_path / "data" / "results" / "dev" / "report.md").read_text()
    assert "Reason atlas (INCONCLUSIVE)" in text and "evaluator: iupac_bases" in text
    config = tmp_path / "bench.json"
    case = json.loads((tmp_path / "data" / "dev" / "c1" / "case.json").read_text())
    gate = case["design"]["depth"] + 1  # every fixture case is now below the gate
    config.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "atlas": {"depth_basis": "design", "min_resolvable_depth": gate},
            }
        )
    )
    cli.main(["--bench-config", str(config), "report", "--split", "dev", "--out-root", data])
    report = json.loads((tmp_path / "data" / "results" / "dev" / "report.json").read_text())
    assert report["sets"][STANDARD]["atlas"]["ladder"]["split_summary"][-1]["expected"] == 2


def test_report_sections_per_set_headline_first_and_rule_on_headline_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    fixture = (("c1", "dupC", "clean"), ("c2", None, STANDARD), ("c3", None, None))

    def fake_run(args: Any) -> tuple[dict[str, Any], int]:
        exact = 1 if args.result_root.name == "hybrid" else 0
        return {"samples": [_sample(n, "normal", "INCONCLUSIVE", exact) for n, e, _ in fixture]}, 0

    monkeypatch.setattr(cli, "evaluate_run", fake_run)
    _split(tmp_path, fixture=fixture)
    for engine in ("ladder", "hybrid"):
        (tmp_path / "data" / "results" / "dev" / engine).mkdir(parents=True)
    data = str(tmp_path / "data")
    cli.main(["evaluate", "--split", "dev", "--engines", "ladder,hybrid", "--out-root", data])
    argv = ["report", "--split", "dev", "--candidate", "hybrid", "--out-root", data]
    assert cli.main(argv) == 0
    report = json.loads((tmp_path / "data" / "results" / "dev" / "report.json").read_text())
    legacy = DEFAULT_BENCH_CONFIG.sets.legacy
    assert report["set_order"] == [STANDARD, "clean", legacy]
    assert {k: v["n_cases"] for k, v in report["sets"].items()} == {
        STANDARD: 1,
        "clean": 1,
        legacy: 1,
    }
    assert report["sets"][STANDARD]["headline"] is True
    assert report["sets"]["clean"]["headline"] is False
    decision = report["decision"]
    assert decision["bench_set"] == STANDARD
    assert decision["profiles"]["ont_amplicon_r10"]["n"] == 1
    text = (tmp_path / "data" / "results" / "dev" / "report.md").read_text()
    first, second = text.index(f"Set `{STANDARD}`"), text.index("Set `clean`")
    assert first < second < text.index(f"Set `{legacy}`")
    assert "headline" in text[first:second]


def test_report_without_headline_cases_skips_the_rule(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _fake_evaluate(cli, monkeypatch, [])
    _split(tmp_path, fixture=(("c1", "dupC", None), ("c2", None, None)))
    for engine in ("ladder", "hybrid"):
        (tmp_path / "data" / "results" / "dev" / engine).mkdir(parents=True)
    data = str(tmp_path / "data")
    cli.main(["evaluate", "--split", "dev", "--engines", "ladder,hybrid", "--out-root", data])
    argv = ["report", "--split", "dev", "--candidate", "hybrid", "--out-root", data]
    assert cli.main(argv) == 0
    report = json.loads((tmp_path / "data" / "results" / "dev" / "report.json").read_text())
    assert report["decision"] is None
    text = (tmp_path / "data" / "results" / "dev" / "report.md").read_text()
    assert f"no `{STANDARD}` cases" in text and "ADOPT" not in text


def test_realism_groups_by_set(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _split(tmp_path, fixture=(("c1", "dupC", "clean"), ("c2", None, STANDARD)))
    monkeypatch.setattr(cli, "case_metrics", lambda *a, **k: {})
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases, cfg: {"n_cases": len(cases)})
    monkeypatch.setattr(cli, "realism_compare", lambda m, t, p, cfg: {})
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    assert cli.main(["realism", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 0
    out = json.loads((tmp_path / "data" / "dev" / "realism.json").read_text())
    assert list(out["sets"]) == [STANDARD, "clean"]
    assert out["sets"]["clean"]["profiles"]["ont_amplicon_r10"]["n_cases"] == 1


def test_only_the_unlocking_rule_can_report_on_test(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Review I1: a second rule registered before unsealing cannot publish a test verdict."""
    cli = _cli(tmp_path, monkeypatch)
    _fake_evaluate(cli, monkeypatch, [])
    _split(tmp_path, "test")
    (tmp_path / "data" / "results" / "test" / "ladder").mkdir(parents=True)
    other = tmp_path / "other.json"
    margin = DEFAULT_BENCH_CONFIG.report.ni_margin * 2
    other.write_text(json.dumps({"schema_version": 1, "report": {"ni_margin": margin}}))
    root = ["--out-root", str(tmp_path / "data")]
    cli.main(["preregister", *root])
    cli.main(["--bench-config", str(other), "preregister", *root])
    assert cli.main(["evaluate", "--split", "test", *root]) == 0
    marker = json.loads((tmp_path / "data" / "test" / "first_evaluation.json").read_text())
    assert marker["rule_sha256"] == rule_sha256(RULE_TEXT)
    assert cli.main(["report", "--split", "test", *root]) == 0
    for action in (["report"], ["evaluate"]):
        with pytest.raises(SystemExit, match="unsealed under rule"):
            cli.main(["--bench-config", str(other), *action, "--split", "test", *root])


def test_realism_index_error_is_a_case_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A case whose metrics raise IndexError is recorded; the split still completes."""
    cli = _cli(tmp_path, monkeypatch)
    _split(tmp_path)

    def metrics(case_dir: Path, *_: Any) -> dict[str, Any]:
        if case_dir.name == "c2":
            raise IndexError("list index out of range")
        return {}

    monkeypatch.setattr(cli, "case_metrics", metrics)
    monkeypatch.setattr(cli, "realism_aggregate", lambda cases, cfg: {"n_cases": len(cases)})
    monkeypatch.setattr(cli, "realism_compare", lambda m, t, p, cfg: {})
    monkeypatch.setattr(cli, "load_targets", lambda: {})
    assert cli.main(["realism", "--split", "dev", "--out-root", str(tmp_path / "data")]) == 0
    out = json.loads((tmp_path / "data" / "dev" / "realism.json").read_text())
    assert out["failures"] == [{"design_id": "c2", "error": "IndexError: list index out of range"}]


def test_report_computes_targets_per_engine_without_a_decision(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Review I4: a per-set root (here `clean` only) still gets its target table."""
    cli = _cli(tmp_path, monkeypatch)
    _fake_evaluate(cli, monkeypatch, [])
    _split(tmp_path, fixture=(("c1", "dupC", "clean"), ("c2", None, "clean")))
    for engine in ("ladder", "hybrid"):
        (tmp_path / "data" / "results" / "dev" / engine).mkdir(parents=True)
    data = str(tmp_path / "data")
    cli.main(["evaluate", "--split", "dev", "--engines", "ladder,hybrid", "--out-root", data])
    assert cli.main(["report", "--split", "dev", "--candidate", "hybrid", "--out-root", data]) == 0
    report = json.loads((tmp_path / "data" / "results" / "dev" / "report.json").read_text())
    assert report["decision"] is None
    hybrid = report["targets"]["hybrid"]
    assert hybrid["clean"]["present"] is True and hybrid["clean"]["pass"] is True
    assert hybrid[STANDARD]["present"] is False and hybrid[STANDARD]["pass"] is None
    assert set(report["targets"]) == {"ladder", "hybrid"}
    text = (tmp_path / "data" / "results" / "dev" / "report.md").read_text()
    assert "descriptive; no decision" in text and "| clean | pooled |" in text
    assert f"{STANDARD}=not present" in text and f"{STANDARD}=False" not in text


def test_report_engines_scored_on_different_cases_exits_cleanly(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _split(tmp_path)
    results = tmp_path / "data" / "results" / "dev"
    samples = {
        "ladder": [_sample("c1", "pathogenic", "PATHOGENIC", 1)],
        "hybrid": [_sample("c2", "normal", "NO_PATHOGENIC_VARIANT_DETECTED", 1)],
    }
    for engine, rows in samples.items():
        (results / engine).mkdir(parents=True)
        (results / engine / "evaluation.json").write_text(json.dumps({"samples": rows}))
    argv = [
        "report",
        "--split",
        "dev",
        "--candidate",
        "hybrid",
        "--out-root",
        str(tmp_path / "data"),
    ]
    with pytest.raises(SystemExit, match="scored on different cases"):
        _cli(tmp_path, monkeypatch).main(argv)


def test_report_records_harness_caller_and_simulator_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from muc_one_span.version import __version__

    cli = _cli(tmp_path, monkeypatch)
    _fake_evaluate(cli, monkeypatch, [])
    split_dir = _split(tmp_path)
    rows = [json.loads(x) for x in (split_dir / "manifest.jsonl").read_text().splitlines()]
    (split_dir / "manifest.jsonl").write_text(
        "".join(json.dumps(r | {"muconeup_version": "0.46.0"}) + "\n" for r in rows)
    )
    engine_dir = tmp_path / "data" / "results" / "dev" / "ladder"
    engine_dir.mkdir(parents=True)
    caller = {"engine": "ladder", "caller_version": "9.9.9", "caller_commit": "abc"}
    (engine_dir / "caller.json").write_text(json.dumps(caller))
    data = str(tmp_path / "data")
    cli.main(["evaluate", "--split", "dev", "--out-root", data])
    assert cli.main(["report", "--split", "dev", "--out-root", data]) == 0
    report = json.loads((engine_dir.parent / "report.json").read_text())
    prov = report["provenance"]
    assert prov["harness_version"] == __version__ and prov["harness_commit"]
    assert prov["callers"] == {"ladder": caller}
    assert prov["muconeup_versions"] == {"0.46.0": 2}
