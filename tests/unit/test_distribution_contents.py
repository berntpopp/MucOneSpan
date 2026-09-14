"""Distribution archives must contain maintained release inputs only."""

from __future__ import annotations

import io
import stat
import tarfile
import zipfile
from pathlib import Path

import pytest

from scripts.check_distribution import validate_sdist, validate_wheel


def write_tar(path: Path, names: list[str], *, symlink: tuple[str, str] | None = None) -> None:
    with tarfile.open(path, "w:gz") as archive:
        for name in names:
            data = b"fixture\n"
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
        if symlink:
            member = tarfile.TarInfo(symlink[0])
            member.type = tarfile.SYMTYPE
            member.linkname = symlink[1]
            archive.addfile(member)


@pytest.mark.parametrize(
    "leak",
    [
        "muc_one_span-0.11.0/tests/results/run/allele_reads.bam",
        "muc_one_span-0.11.0/tests/data/generated/sample/reads.fastq",
        "muc_one_span-0.11.0/src/muc_one_span/reads.fastq.gz",
        "muc_one_span-0.11.0/src/muc_one_span/data/reference/reference_ladder.fa.fai",
        "muc_one_span-0.11.0/src/muc_one_span/__pycache__/cli.cpython-312.pyc",
        "muc_one_span-0.11.0/.planning/reviews/private-response.json",
    ],
)
def test_sdist_rejects_generated_and_private_content(tmp_path: Path, leak: str) -> None:
    archive = tmp_path / "package.tar.gz"
    write_tar(
        archive,
        [
            "muc_one_span-0.11.0/pyproject.toml",
            "muc_one_span-0.11.0/src/muc_one_span/version.py",
            leak,
        ],
    )
    with pytest.raises(ValueError, match=leak):
        validate_sdist(archive)


def test_sdist_accepts_explicit_maintained_source_inventory(tmp_path: Path) -> None:
    archive = tmp_path / "package.tar.gz"
    write_tar(
        archive,
        [
            "muc_one_span-0.11.0/pyproject.toml",
            "muc_one_span-0.11.0/README.md",
            "muc_one_span-0.11.0/LICENSE",
            "muc_one_span-0.11.0/src/muc_one_span/version.py",
            "muc_one_span-0.11.0/tests/unit/test_cli.py",
            "muc_one_span-0.11.0/tests/integration/test_pipeline.py",
            "muc_one_span-0.11.0/scripts/check_distribution.py",
            "muc_one_span-0.11.0/docs/index.md",
            "muc_one_span-0.11.0/examples/runtime-settings.json",
        ],
    )
    validate_sdist(archive)


def test_sdist_rejects_symlink_and_parent_traversal(tmp_path: Path) -> None:
    symlink_archive = tmp_path / "symlink.tar.gz"
    write_tar(
        symlink_archive,
        ["muc_one_span-0.11.0/pyproject.toml"],
        symlink=("muc_one_span-0.11.0/src/muc_one_span/data-link", "../../private"),
    )
    with pytest.raises(ValueError, match="symbolic or hard link"):
        validate_sdist(symlink_archive)

    traversal_archive = tmp_path / "traversal.tar.gz"
    write_tar(traversal_archive, ["muc_one_span-0.11.0/../secret"])
    with pytest.raises(ValueError, match="unsafe archive path"):
        validate_sdist(traversal_archive)


def test_wheel_rejects_generated_reads_and_symlinks(tmp_path: Path) -> None:
    wheel = tmp_path / "package.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("muc_one_span/version.py", '__version__ = "0.11.0"\n')
        archive.writestr("muc_one_span/tests/data/generated/reads.fastq", "@r\nA\n+\nI\n")
    with pytest.raises(ValueError, match=r"reads\.fastq"):
        validate_wheel(wheel)

    symlink_wheel = tmp_path / "symlink.whl"
    with zipfile.ZipFile(symlink_wheel, "w") as archive:
        info = zipfile.ZipInfo("muc_one_span/data-link")
        info.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(info, "../../private")
    with pytest.raises(ValueError, match="symbolic link"):
        validate_wheel(symlink_wheel)
