"""benchsim.provenance: harness and caller commits, caller.json round trip."""

from pathlib import Path
from typing import Any

from muc_one_span.benchsim.provenance import (
    caller_record,
    git_commit,
    read_caller,
    report_provenance,
    write_caller,
)
from muc_one_span.version import __version__


def test_git_commit_uses_the_runner_and_tolerates_no_git(tmp_path: Path) -> None:
    seen: list[Any] = []

    def ok(args: list[str], cwd: str) -> str:
        seen.append((args, cwd))
        return "abc123\n"

    def missing(args: list[str], cwd: str) -> str:
        raise FileNotFoundError("git")

    def failed(args: list[str], cwd: str) -> str:
        raise RuntimeError("not a repository")

    assert git_commit(tmp_path, ok) == "abc123"
    assert seen == [(["git", "rev-parse", "HEAD"], str(tmp_path))]
    assert git_commit(tmp_path, missing) is None and git_commit(tmp_path, failed) is None
    assert git_commit(tmp_path, lambda args, cwd: "\n") is None


def test_caller_json_round_trip_and_report_summary(tmp_path: Path) -> None:
    engine_dir = tmp_path / "ladder"
    record = caller_record("ladder", "abc")
    assert record == {"engine": "ladder", "caller_version": __version__, "caller_commit": "abc"}
    write_caller(engine_dir, record)
    assert read_caller(engine_dir) == record
    assert read_caller(tmp_path / "hybrid") is None
    cases = {"a": {"muconeup_version": "0.46.0"}, "b": {"muconeup_version": "0.45.0"}, "c": {}}
    prov = report_provenance("def", {"ladder": engine_dir}, cases)
    assert prov["harness_commit"] == "def" and prov["callers"] == {"ladder": record}
    assert prov["muconeup_versions"] == {"0.45.0": 1, "0.46.0": 1, "None": 1}
