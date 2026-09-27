"""Registering decision rule v6 on a split already pre-registered under v5 (task 15o).

Mirrors a sealed `test` split generated and pre-registered before target scopes
existed: its bench-config file names every target without a ``scope``, and its ledger
holds one earlier rule line. v6 is appended with a derived config (``scope: "pooled"``
on the INCONCLUSIVE targets), then must be accepted by ``evaluate`` and ``report``;
nothing a scope touches may change the generation hash.
"""

import copy
import json
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, load_bench_config
from muc_one_span.benchsim.preregistration import preregister
from muc_one_span.benchsim.report import RULE_TEXT, rule_sha256
from tests.unit.benchsim.test_benchsim_cli import _cli
from tests.unit.benchsim.test_benchsim_report_cli import _fake_evaluate, _split

PRE_V6_RULE = "MucSim-Bench decision rule v5 (stand-in for the earlier registered text)"


def _pre_v6_config() -> dict[str, Any]:
    """The effective settings as a pre-v6 file wrote them: targets without a scope."""
    data = DEFAULT_BENCH_CONFIG.to_dict()
    data["design"]["split_sizes"]["test"] = 150  # a non-default split size, as on test150b
    for metrics in data["targets"]["by_set"].values():
        for target in metrics.values():
            target.pop("scope")
    return data


def _v6_config(pre_v6: dict[str, Any]) -> dict[str, Any]:
    """The documented migration: the same file with INCONCLUSIVE scoped to the pooled set."""
    data = copy.deepcopy(pre_v6)
    for metrics in data["targets"]["by_set"].values():
        metrics["inconclusive_rate"]["scope"] = "pooled"
    return data


def _write(path: Path, data: dict[str, Any]) -> Path:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    return path


def test_v6_registers_next_to_v5_and_is_accepted_by_evaluate_and_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli = _cli(tmp_path, monkeypatch)
    _fake_evaluate(cli, monkeypatch, [])
    _split(tmp_path, "test")
    for engine in ("ladder", "hybrid"):
        (tmp_path / "data" / "results" / "test" / engine).mkdir(parents=True)
    ledger = tmp_path / "data" / "test" / "preregistration.jsonl"
    preregister(PRE_V6_RULE, ledger)
    first_line = ledger.read_text()
    old = _write(tmp_path / "bench-pre-v6.json", _pre_v6_config())
    new = _write(tmp_path / "bench-v6.json", _v6_config(_pre_v6_config()))
    old_cfg, new_cfg = load_bench_config(old), load_bench_config(new)

    # Scopes are not generation-shaping: every case stays reusable and verifiable.
    for name, definition in new_cfg.sets.definitions.items():
        for profile in definition.profiles:
            assert new_cfg.generation_sha256(name, profile) == old_cfg.generation_sha256(
                name, profile
            )
    # The derived file yields exactly the default v6 rule (same report/headline settings).
    root = ["--out-root", str(tmp_path / "data")]
    assert cli.main(["--bench-config", str(new), "preregister", *root]) == 0
    lines = ledger.read_text().splitlines()
    assert ledger.read_text().startswith(first_line) and len(lines) == 2
    assert json.loads(lines[1])["sha256"] == rule_sha256(RULE_TEXT)

    # A scope-less (pre-v6) file keeps the old per-profile semantics: another rule.
    with pytest.raises(SystemExit, match="not pre-registered"):
        cli.main(["--bench-config", str(old), "report", "--split", "test", *root])

    evaluate = ["--bench-config", str(new), "evaluate", "--engines", "ladder,hybrid"]
    assert cli.main([*evaluate, "--split", "test", *root]) == 0
    marker = json.loads((tmp_path / "data" / "test" / "first_evaluation.json").read_text())
    assert marker["rule_sha256"] == rule_sha256(RULE_TEXT)
    report_cmd = ["--bench-config", str(new), "report", "--candidate", "hybrid"]
    assert cli.main([*report_cmd, "--split", "test", *root]) == 0
    report = json.loads((tmp_path / "data" / "results" / "test" / "report.json").read_text())
    assert report["decision"]["rule_sha256"] == rule_sha256(RULE_TEXT)
    assert report["preregistration"]["sha256"] == rule_sha256(RULE_TEXT)
    assert report["bench_config_sha256"] == new_cfg.sha256()
    markdown = (tmp_path / "data" / "results" / "test" / "report.md").read_text()
    assert f"Rule sha256: `{rule_sha256(RULE_TEXT)}`" in markdown
    # Once unsealed under v6, the earlier rule can no longer publish a verdict.
    with pytest.raises(SystemExit):
        cli.main(["--bench-config", str(old), "report", "--split", "test", *root])
