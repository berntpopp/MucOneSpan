"""Test unsealing is bound to the one pre-registered rule that unlocked it (review I1)."""

import json
from pathlib import Path

import pytest

from muc_one_span.benchsim.preregistration import (
    FIRST_EVALUATION,
    mark_first_evaluation,
    preregister,
    require_preregistered,
    rule_sha256,
    unlocking_rule_sha256,
)


def _ledger(tmp_path: Path) -> Path:
    path = tmp_path / "test" / "preregistration.jsonl"
    preregister("rule A", path)
    preregister("rule B", path)
    return path


def test_first_evaluation_records_the_unlocking_rule(tmp_path: Path) -> None:
    path = _ledger(tmp_path)
    assert unlocking_rule_sha256(path) is None
    mark_first_evaluation(path, rule_sha256("rule A"))
    marker = json.loads((path.parent / FIRST_EVALUATION).read_text())
    assert marker["rule_sha256"] == rule_sha256("rule A") == unlocking_rule_sha256(path)


def test_only_the_unlocking_rule_is_accepted_after_unsealing(tmp_path: Path) -> None:
    path = _ledger(tmp_path)
    require_preregistered(path, "rule B")  # both registered, test still sealed
    mark_first_evaluation(path, rule_sha256("rule A"))
    assert require_preregistered(path, "rule A")["sha256"] == rule_sha256("rule A")
    with pytest.raises(PermissionError, match="unsealed under rule"):
        require_preregistered(path, "rule B")


def test_marking_again_under_another_rule_is_refused(tmp_path: Path) -> None:
    path = _ledger(tmp_path)
    stamp = mark_first_evaluation(path, rule_sha256("rule A"))
    assert mark_first_evaluation(path, rule_sha256("rule A")) == stamp
    with pytest.raises(PermissionError, match="unsealed under rule"):
        mark_first_evaluation(path, rule_sha256("rule B"))


def test_a_marker_without_a_rule_hash_fails_closed(tmp_path: Path) -> None:
    path = _ledger(tmp_path)
    (path.parent / FIRST_EVALUATION).write_text(json.dumps({"evaluated_at": "then"}) + "\n")
    with pytest.raises(PermissionError, match="records no unlocking rule"):
        require_preregistered(path, "rule A")
