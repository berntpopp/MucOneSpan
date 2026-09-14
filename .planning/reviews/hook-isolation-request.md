Review this bounded test-only release fix with Claude Fable 5.1. No production behavior changed after final release review. Base reviewed release commit 5f62cf296a4036c9dd86200438684b42e5d2e87d. Push hook ran all 694 unit tests successfully but failed because the existing file-size tests inherited Git hook environment and staged workflow.yml in the real invoking index. That accidental index entry was removed. New two-context regression runs the existing tests under isolated fake invoking-repository GIT_INDEX_FILE or GIT_DIR/GIT_WORK_TREE and verifies byte-identical index afterward. Both new tests failed before the fix; five focused tests pass after adding an autouse fixture scoped to test_file_size.py that unsets names returned by git rev-parse --local-env-vars. Check correctness, isolation, and whether any gate is weakened. Read-only review; concise blocking findings or no blockers.

Diff:
diff --git a/tests/unit/test_file_size.py b/tests/unit/test_file_size.py
index 79d667e..50dc056 100644
--- a/tests/unit/test_file_size.py
+++ b/tests/unit/test_file_size.py
@@ -9,6 +9,14 @@ import pytest
 CHECKER = Path(__file__).resolve().parents[2] / "scripts/check_file_size.py"


+@pytest.fixture(autouse=True)
+def isolate_temporary_git_repository(monkeypatch: pytest.MonkeyPatch) -> None:
+    """Prevent hook-local Git state from redirecting temporary repository commands."""
+    names = subprocess.check_output(["git", "rev-parse", "--local-env-vars"], text=True)
+    for name in names.splitlines():
+        monkeypatch.delenv(name, raising=False)
+
+
 @pytest.mark.parametrize("lines, expected", [(649, 0), (650, 1)])
 def test_strict_boundary_includes_untracked_code(tmp_path: Path, lines: int, expected: int) -> None:
     subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)

New test:
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
