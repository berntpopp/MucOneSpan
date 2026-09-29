"""Registering decision rule v7 on a sealed split already pre-registered under v6 (task 15n).

Mirrors test150b: a bench-config file derived for v6 (INCONCLUSIVE targets scoped to the
pooled set, the clean ceiling still 0.10) and a ledger holding the earlier rule lines. The
documented v7 migration changes only ``targets.by_set.clean.inconclusive_rate.threshold``
to 0.15; the derived file must give exactly the default v7 rule, keep every generation
hash, and be accepted by ``evaluate`` and ``report``.
"""

import copy
import json
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, load_bench_config
from muc_one_span.benchsim.preregistration import preregister
from muc_one_span.benchsim.report import RULE_TEXT, rule_sha256, rule_text
from tests.unit.benchsim.test_benchsim_cli import _cli
from tests.unit.benchsim.test_benchsim_report_cli import _fake_evaluate, _split
from tests.unit.benchsim.test_rule_v6_migration import _write

V6_CLEAN_INCONCLUSIVE = 0.10  # the v6 default the owner relaxed (decision 2026-09-28)


def _v6_config() -> dict[str, Any]:
    """The effective settings as a v6-derived file holds them (clean ceiling 0.10)."""
    data = DEFAULT_BENCH_CONFIG.to_dict()
    data["design"]["split_sizes"]["test"] = 150  # a non-default split size, as on test150b
    data["targets"]["by_set"]["clean"]["inconclusive_rate"]["threshold"] = V6_CLEAN_INCONCLUSIVE
    return data


def _v7_config(v6: dict[str, Any]) -> dict[str, Any]:
    """The documented migration: only the clean INCONCLUSIVE threshold changes."""
    data = copy.deepcopy(v6)
    default = DEFAULT_BENCH_CONFIG.targets.by_set["clean"]["inconclusive_rate"].threshold
    data["targets"]["by_set"]["clean"]["inconclusive_rate"]["threshold"] = default
    return data


def test_v7_differs_from_a_v6_file_only_by_the_clean_inconclusive_ceiling(
    tmp_path: Path,
) -> None:
    v6 = load_bench_config(_write(tmp_path / "bench-v6.json", _v6_config()))
    text = rule_text(v6.report, v6.sets.headline, v6.targets)
    assert text != RULE_TEXT
    v6_clause = "inconclusive_rate <= 0.1 (pooled only"
    assert text.count(v6_clause) == 1
    assert text.replace(v6_clause, "inconclusive_rate <= 0.15 (pooled only") == RULE_TEXT


def test_v7_registers_next_to_v6_and_is_accepted_by_evaluate_and_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _fake_evaluate(cli, monkeypatch, [])
    _split(tmp_path, "test")
    for engine in ("ladder", "hybrid"):
        (tmp_path / "data" / "results" / "test" / engine).mkdir(parents=True)
    ledger = tmp_path / "data" / "test" / "preregistration.jsonl"
    old = _write(tmp_path / "bench-v6.json", _v6_config())
    new = _write(tmp_path / "bench-v7.json", _v7_config(_v6_config()))
    old_cfg, new_cfg = load_bench_config(old), load_bench_config(new)
    preregister(rule_text(old_cfg.report, old_cfg.sets.headline, old_cfg.targets), ledger)
    first_line = ledger.read_text()

    for name, definition in new_cfg.sets.definitions.items():
        for profile in definition.profiles:
            assert new_cfg.generation_sha256(name, profile) == old_cfg.generation_sha256(
                name, profile
            )
    root = ["--out-root", str(tmp_path / "data")]
    assert cli.main(["--bench-config", str(new), "preregister", *root]) == 0
    lines = ledger.read_text().splitlines()
    assert ledger.read_text().startswith(first_line) and len(lines) == 2
    assert json.loads(lines[1])["sha256"] == rule_sha256(RULE_TEXT)
    assert json.loads(lines[0])["sha256"] != rule_sha256(RULE_TEXT)

    evaluate = ["--bench-config", str(new), "evaluate", "--engines", "ladder,hybrid"]
    assert cli.main([*evaluate, "--split", "test", *root]) == 0
    marker = json.loads((tmp_path / "data" / "test" / "first_evaluation.json").read_text())
    assert marker["rule_sha256"] == rule_sha256(RULE_TEXT)
    report_cmd = ["--bench-config", str(new), "report", "--candidate", "hybrid"]
    assert cli.main([*report_cmd, "--split", "test", *root]) == 0
    report = json.loads((tmp_path / "data" / "results" / "test" / "report.json").read_text())
    assert report["decision"]["rule_sha256"] == rule_sha256(RULE_TEXT)
    # Once unsealed under v7, the v6 file can no longer publish a verdict.
    with pytest.raises(SystemExit):
        cli.main(["--bench-config", str(old), "report", "--split", "test", *root])
