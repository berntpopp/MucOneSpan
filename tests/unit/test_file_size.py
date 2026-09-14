"""Exercise the repository size gate against real temporary Git worktrees."""

import subprocess
import sys
from pathlib import Path

import pytest

CHECKER = Path(__file__).resolve().parents[2] / "scripts/check_file_size.py"


@pytest.mark.parametrize("lines, expected", [(649, 0), (650, 1)])
def test_strict_boundary_includes_untracked_code(tmp_path: Path, lines: int, expected: int) -> None:
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    source = tmp_path / "new module.py"
    source.write_text("# line\n" * (lines - 1) + "# final line", encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(CHECKER)], cwd=tmp_path, text=True, capture_output=True
    )
    assert result.returncode == expected, result.stderr
    if expected:
        assert f"new module.py: {lines}" in result.stdout


def test_ignores_generated_files_but_checks_tracked_code(tmp_path: Path) -> None:
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / ".gitignore").write_text("cache/\n", encoding="utf-8")
    (tmp_path / "cache").mkdir()
    (tmp_path / "cache/generated.py").write_text("x\n" * 900, encoding="utf-8")
    (tmp_path / "uv.lock").write_text("x\n" * 900, encoding="utf-8")
    workflow = tmp_path / "workflow.yml"
    workflow.write_text("# line\n" * 650, encoding="utf-8")
    subprocess.run(["git", "add", "workflow.yml"], cwd=tmp_path, check=True)
    result = subprocess.run(
        [sys.executable, str(CHECKER)], cwd=tmp_path, text=True, capture_output=True
    )
    assert result.returncode == 1
    assert "workflow.yml: 650" in result.stdout
    assert "cache" not in result.stdout
    assert "uv.lock" not in result.stdout
