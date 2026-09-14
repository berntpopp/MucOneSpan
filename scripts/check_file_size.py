#!/usr/bin/env python3
"""Enforce fewer than 650 physical lines in authored code and configuration.

Scan tracked and new nonignored files, including files outside the Python package.
Lockfiles, prose and scientific data are not source code and are excluded by type.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

MAX_LINES = 649
SOURCE_SUFFIXES = {".py", ".sh", ".yml", ".yaml", ".toml", ".css", ".js", ".ts", ".html", ".j2"}
SOURCE_NAMES = {"Makefile", "Dockerfile"}


def main() -> int:
    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        check=True,
        capture_output=True,
    )
    failures = []
    checked = 0
    for name in sorted(set(result.stdout.decode().split("\0")) - {""}):
        path = Path(name)
        if path.suffix not in SOURCE_SUFFIXES and path.name not in SOURCE_NAMES:
            continue
        if not path.is_file():  # Tracked files may have been deleted in the worktree.
            continue
        lines = len(path.read_bytes().splitlines())
        checked += 1
        if lines > MAX_LINES:
            failures.append(f"{name}: {lines} lines (maximum {MAX_LINES}); split by responsibility")
    if failures:
        print("\n".join(failures))
        return 1
    print(
        f"File size check passed: {checked} source/configuration files, maximum {MAX_LINES} lines"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
