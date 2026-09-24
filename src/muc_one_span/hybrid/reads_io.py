"""Streaming read input for the hybrid engine and the versions of its optional packages.

FASTQ and gzipped FASTQ are parsed in Python, so the hybrid path runs no external tool
for FASTQ input. A BAM is streamed through ``samtools fastq`` (argument list, no shell)
via ``tools.run_tool_iter``, keeping primary records only.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from importlib import metadata
from pathlib import Path

from muc_one_span.hybrid.spans import ReadRecord

# SAM FLAG mask excluded by ``samtools fastq -F``: secondary (0x100) | supplementary
# (0x800) alignments, so each molecule is read once (format definition, not a tunable).
NON_PRIMARY_FLAGS = "0x900"
FASTQ_HEADER = "@"
# Distribution name of each POA backend's Python package (for version provenance).
BACKEND_PACKAGES = {"pyabpoa": "pyabpoa", "pyspoa": "pyspoa"}
UNKNOWN_VERSION = "unknown"


def parse_fastq(lines: Iterable[str]) -> Iterator[ReadRecord]:
    """Yield four-line FASTQ records; a truncated or malformed record is an error."""
    it = iter(lines)
    for header in it:
        if not header.strip():
            continue
        if not header.startswith(FASTQ_HEADER):
            raise ValueError(f"malformed FASTQ header: {header.strip()[:40]!r}")
        try:
            seq, _plus, qual = next(it), next(it), next(it)
        except StopIteration:
            raise ValueError("truncated FASTQ record") from None
        name = header[len(FASTQ_HEADER) :].split()
        yield ReadRecord(name[0] if name else "", seq.strip().upper(), qual.strip())


def read_input(path: Path) -> Iterator[ReadRecord]:
    """Stream FASTQ(.gz) records, or the primary reads of a BAM via ``samtools fastq``."""
    if path.suffix == ".bam":
        from muc_one_span.tools import run_tool_iter

        yield from parse_fastq(
            run_tool_iter(["samtools", "fastq", "-F", NON_PRIMARY_FLAGS, str(path)])
        )
    elif path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            yield from parse_fastq(handle)
    else:
        with path.open() as handle:
            yield from parse_fastq(handle)


def extra_versions(backend: str) -> dict[str, str]:
    """Versions of the optional packages that produced the consensus."""
    out = {}
    for pkg in ("edlib", BACKEND_PACKAGES[backend]):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = UNKNOWN_VERSION
    return out
