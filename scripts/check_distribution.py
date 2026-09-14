#!/usr/bin/env python3
"""Validate release archive contents and smoke-test the wheel."""

from __future__ import annotations

import stat
import subprocess
import tarfile
import zipfile
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

SDIST_DIRECTORIES = ("src/", "tests/unit/", "tests/integration/", "scripts/", "docs/", "examples/")
SDIST_FILES = {
    ".gitignore",
    ".python-version",
    "CHANGELOG.md",
    "CITATION.cff",
    "LICENSE",
    "Makefile",
    "PKG-INFO",
    "README.md",
    "mkdocs.yml",
    "pyproject.toml",
    "tests/__init__.py",
    "tests/conftest.py",
    "uv.lock",
}
FORBIDDEN_PARTS = {".planning", "__pycache__", "htmlcov", "results", "site"}
FORBIDDEN_SUFFIXES = {
    ".amb",
    ".ann",
    ".bai",
    ".bam",
    ".bwt",
    ".crai",
    ".cram",
    ".fai",
    ".fastq",
    ".fq",
    ".mmi",
    ".pac",
    ".pyc",
    ".pyo",
    ".sa",
}
FORBIDDEN_COMPOUND_SUFFIXES = (".fastq.gz", ".fq.gz")
REQUIRED_WHEEL_FILES = {
    "muc_one_span/version.py",
    "muc_one_span/data/reference/reference_ladder.fa",
    "muc_one_span/data/repeats/repeats.json",
    "muc_one_span/templates/report.html.j2",
}


def _unsafe_path(name: str) -> bool:
    path = PurePosixPath(name)
    return not name or "\\" in name or path.is_absolute() or ".." in path.parts


def _forbidden_path(path: PurePosixPath) -> bool:
    name = path.name.lower()
    return (
        bool(set(path.parts) & FORBIDDEN_PARTS)
        or path.suffix.lower() in FORBIDDEN_SUFFIXES
        or name.endswith(FORBIDDEN_COMPOUND_SUFFIXES)
    )


def _raise_violations(kind: str, violations: list[str]) -> None:
    if violations:
        shown = "\n".join(f"- {value}" for value in violations[:25])
        remaining = f"\n- ... {len(violations) - 25} more" if len(violations) > 25 else ""
        raise ValueError(
            f"Unsafe or unexpected {kind} contents ({len(violations)}):\n{shown}{remaining}"
        )


def validate_sdist(path: Path) -> None:
    """Reject unsafe links, generated data, and files outside the sdist allowlist."""
    violations: list[str] = []
    seen: set[str] = set()
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        roots = {PurePosixPath(member.name).parts[0] for member in members if member.name}
        if len(roots) != 1:
            violations.append(f"sdist must have one archive root, found {sorted(roots)}")
        root = next(iter(roots), "")
        for member in members:
            if _unsafe_path(member.name):
                violations.append(f"unsafe archive path: {member.name}")
                continue
            if member.issym() or member.islnk():
                violations.append(f"symbolic or hard link: {member.name} -> {member.linkname}")
                continue
            if not (member.isfile() or member.isdir()):
                violations.append(f"unsupported archive member type: {member.name}")
                continue
            relative = PurePosixPath(*PurePosixPath(member.name).parts[1:])
            if member.isdir() or not relative.parts:
                continue
            value = relative.as_posix()
            seen.add(value)
            allowed = value in SDIST_FILES or value.startswith(SDIST_DIRECTORIES)
            if not allowed or _forbidden_path(relative):
                violations.append(member.name)
        for required in ("pyproject.toml", "README.md", "src/muc_one_span/version.py"):
            if required not in seen:
                violations.append(f"missing required file below {root}: {required}")
    _raise_violations("sdist", violations)


def validate_wheel(path: Path) -> None:
    """Reject unsafe/generated wheel members and require runtime resources."""
    violations: list[str] = []
    seen: set[str] = set()
    with zipfile.ZipFile(path) as archive:
        for member in archive.infolist():
            if _unsafe_path(member.filename):
                violations.append(f"unsafe archive path: {member.filename}")
                continue
            mode = member.external_attr >> 16
            if stat.S_ISLNK(mode):
                violations.append(f"symbolic link: {member.filename}")
                continue
            if member.is_dir():
                continue
            relative = PurePosixPath(member.filename)
            seen.add(relative.as_posix())
            top = relative.parts[0] if relative.parts else ""
            allowed = top == "muc_one_span" or top.endswith(".dist-info")
            if not allowed or _forbidden_path(relative):
                violations.append(member.filename)
        for required in REQUIRED_WHEEL_FILES - seen:
            violations.append(f"missing required wheel file: {required}")
    _raise_violations("wheel", violations)


def main() -> None:
    wheels = list(Path("dist").glob("*.whl"))
    sdists = list(Path("dist").glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise SystemExit("Expected exactly one wheel and one sdist in dist")
    wheel = wheels[0].resolve()
    sdist = sdists[0].resolve()
    try:
        validate_wheel(wheel)
        validate_sdist(sdist)
    except (OSError, tarfile.TarError, zipfile.BadZipFile, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    with TemporaryDirectory(prefix="muconespan-package-") as directory:
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
                "from importlib.metadata import distribution; "
                "from importlib.resources import files; "
                "from muc_one_span import __version__; "
                "from muc_one_span.config import load_repeat_dictionary; "
                "from muc_one_span.report import generate_report; "
                "from pathlib import Path; "
                "package = distribution('muc_one_span'); "
                "assert package.version == __version__; "
                "assert {(e.name, e.value) for e in package.entry_points "
                "if e.group == 'console_scripts'} == {('muconespan', 'muc_one_span.cli:main')}; "
                "assert load_repeat_dictionary().repeats; "
                "assert files('muc_one_span').joinpath('templates/report.html.j2').is_file(); "
                "assert files('muc_one_span.data.reference').joinpath('reference_ladder.fa').is_file(); "
                "generate_report({'alleles': {}, 'classifications': {}}, Path('report.html')); "
                "assert '--bg: #ffffff' in Path('report.html').read_text()",
            ],
            cwd=root,
            check=True,
        )
        subprocess.run([str(root / "venv/bin/muconespan"), "--help"], cwd=root, check=True)
        subprocess.run([str(root / "venv/bin/muconespan"), "--version"], cwd=root, check=True)
    print("Wheel/sdist contents, wheel installation, CLI and bundled resources passed")


if __name__ == "__main__":
    main()
