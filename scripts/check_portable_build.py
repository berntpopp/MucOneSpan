#!/usr/bin/env python3
"""Fail when an installed compiled extension uses AVX-family (VEX/EVEX) instructions.

``pyabpoa`` is published as a source distribution only, and its ``setup.py``
compiles with ``-march=native`` on x86-64 Linux unless the build environment sets
``SSE4`` (``-msse4.1``) or ``SSE2`` (``-msse2``). A native build crashes with
"Illegal instruction" on any CPU that lacks the build host's instruction set:
cached CI wheels, container images and shared cluster environments are affected.

The project builds pyabpoa with ``SSE4=1`` (``[tool.uv.extra-build-variables]`` in
``pyproject.toml`` and the Docker builder stage). This guard disassembles each
named extension with ``objdump`` (GNU binutils) and rejects AVX, AVX2 and AVX-512
instructions, which the SSE4.1 baseline never emits. It is x86-64 specific: on
aarch64 pyabpoa's ``setup.py`` already targets the portable ``armv8-a+simd``
baseline, so the guard reports a skip and exits 0 on other architectures.

Exit status: 0 portable, 1 non-portable instructions found, 2 guard error.
"""

from __future__ import annotations

import argparse
import importlib.machinery
import importlib.util
import platform
import re
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

# AT&T-syntax objdump format definitions, not tunables. Every VEX/EVEX-encoded
# (AVX-family) mnemonic starts with "v"; AVX-512 mask instructions start with
# "k" and use %k registers; 256/512-bit operands use %ymm/%zmm registers.
INSTRUCTION_LINE = re.compile(r"^\s*[0-9a-f]+:\t(?P<text>\S.*?)\s*$")
AVX_FAMILY = re.compile(r"^(?:v|k[a-z])\S*\s|%[yz]mm\d|%k[0-7]\b")
REPORTED_EXAMPLES = 5
X86_64_MACHINES = frozenset({"x86_64", "amd64", "AMD64"})
REMEDY = (
    "Rebuild portably: `uv cache clean pyabpoa`, then reinstall with SSE4=1 "
    "(uv applies it from pyproject.toml; pip needs `SSE4=1 pip install --no-cache-dir`)."
)


class GuardError(RuntimeError):
    """The guard could not inspect an extension."""


def locate_extension(module: str) -> Path:
    """Return the compiled extension file of an installed top-level module."""
    spec = importlib.util.find_spec(module)
    if spec is None or spec.origin is None:
        raise GuardError(f"{module} is not installed in {sys.executable}")
    origin = Path(spec.origin)
    if not origin.name.endswith(tuple(importlib.machinery.EXTENSION_SUFFIXES)):
        raise GuardError(f"{module} is not a compiled extension: {origin}")
    return origin


def find_non_portable_instructions(listing: str) -> list[str]:
    """Return the AVX-family instructions in an ``objdump -d`` listing."""
    found = []
    for line in listing.splitlines():
        match = INSTRUCTION_LINE.match(line)
        if match and AVX_FAMILY.search(match.group("text") + " "):
            found.append(match.group("text"))
    return found


def disassemble(extension: Path, objdump: str) -> str:
    """Disassemble an extension with objdump, raising GuardError on failure."""
    command = [objdump, "-d", "--no-show-raw-insn", str(extension)]
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise GuardError(f"objdump not found ({objdump}); install GNU binutils") from exc
    except subprocess.CalledProcessError as exc:
        raise GuardError(f"objdump failed on {extension}: {exc.stderr}") from exc
    return result.stdout


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("modules", nargs="+", help="installed extension modules to inspect")
    parser.add_argument("--objdump", default="objdump", help="objdump executable")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Check every named module; return the exit status."""
    args = parse_args(argv)
    if platform.machine() not in X86_64_MACHINES:
        print(f"skipped: the AVX guard is x86-64 specific (machine {platform.machine()})")
        return 0
    status = 0
    for module in args.modules:
        try:
            extension = locate_extension(module)
            offending = find_non_portable_instructions(disassemble(extension, args.objdump))
        except GuardError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        if not offending:
            print(f"{module}: portable (no AVX-family instructions in {extension})")
            continue
        status = 1
        print(f"{module}: {len(offending)} AVX-family instruction(s) in {extension}")
        for text in offending[:REPORTED_EXAMPLES]:
            print(f"  {text}")
        print(REMEDY)
    return status


if __name__ == "__main__":
    sys.exit(main())
