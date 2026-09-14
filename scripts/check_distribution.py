#!/usr/bin/env python3
"""Smoke-test the built wheel outside the checkout in an isolated environment."""

from __future__ import annotations

import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory


def main() -> None:
    wheels = list(Path("dist").glob("*.whl"))
    if len(wheels) != 1:
        raise SystemExit("Expected exactly one wheel in dist; remove stale artifacts first")
    wheel = wheels[0].resolve()
    with TemporaryDirectory(prefix="pacmuci-package-") as directory:
        root = Path(directory)
        python = root / "venv/bin/python"
        subprocess.run(["uv", "venv", str(root / "venv")], check=True)
        subprocess.run(
            ["uv", "pip", "install", "--python", str(python), f"{wheel}[report]"], check=True
        )
        subprocess.run(
            [
                str(python),
                "-c",
                "from importlib.resources import files; "
                "from open_pacmuci.config import load_repeat_dictionary; "
                "from open_pacmuci.report import generate_report; "
                "from pathlib import Path; "
                "assert load_repeat_dictionary().repeats; "
                "assert files('open_pacmuci').joinpath('templates/report.html.j2').is_file(); "
                "assert files('open_pacmuci.data.reference').joinpath('reference_ladder.fa').is_file(); "
                "generate_report({'alleles': {}, 'classifications': {}}, Path('report.html')); "
                "assert '--bg: #ffffff' in Path('report.html').read_text()",
            ],
            cwd=root,
            check=True,
        )
        subprocess.run([str(root / "venv/bin/open-pacmuci"), "--help"], cwd=root, check=True)
    print("Wheel installation, CLI and bundled resources passed")


if __name__ == "__main__":
    main()
