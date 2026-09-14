"""Repository-gate tests must not mutate a Git hook's invoking repository."""

import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("hook_context", ["index", "repository"])
def test_size_tests_preserve_invoking_git_state(tmp_path: Path, hook_context: str) -> None:
    env = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    invoking = tmp_path / "invoking"
    invoking.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=invoking, env=env, check=True)
    (invoking / "keep.txt").write_text("preserve this staged file\n")
    subprocess.run(["git", "add", "keep.txt"], cwd=invoking, env=env, check=True)
    index = invoking / ".git/index"
    before = index.read_bytes()
    if hook_context == "index":
        env["GIT_INDEX_FILE"] = str(index)
    else:
        env.update(GIT_DIR=str(invoking / ".git"), GIT_WORK_TREE=str(invoking))
    target = Path(__file__).with_name("test_file_size.py")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(target), "--no-cov", "-q"],
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert index.read_bytes() == before, "temporary tests changed the invoking Git index"
