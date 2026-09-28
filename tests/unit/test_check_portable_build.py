"""Exercise the portable-build guard with mocked disassembly (no objdump needed)."""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

SCRIPT = Path(__file__).parents[2] / "scripts" / "check_portable_build.py"

# objdump -d --no-show-raw-insn excerpts: an SSE4.1 build and the AVX families
# that a -march=native build emits on newer CPUs.
SSE4_LISTING = """
pyabpoa.so:     file format elf64-x86-64

Disassembly of section .text:

0000000000001000 <abpoa_align>:
    1000:\tpush   %rbp
    1001:\tpminsd %xmm1,%xmm0
    1006:\tpblendvb %xmm0,%xmm2,%xmm3
    100b:\tmovdqa %xmm0,(%rdi)
    100f:\tret
"""
AVX2_LISTING = SSE4_LISTING + "    1010:\tvpmaxsd %ymm1,%ymm2,%ymm0\n"
AVX512_LISTING = SSE4_LISTING + "    1020:\tvpaddd %zmm1,%zmm2,%zmm0{%k1}\n"
MASK_ONLY_LISTING = SSE4_LISTING + "    1030:\tkmovw  %k1,%eax\n"
VEX_XMM_LISTING = SSE4_LISTING + "    1040:\tvmovdqa %xmm0,(%rdi)\n"


def load_script() -> ModuleType:
    spec = importlib.util.spec_from_file_location("check_portable_build_script", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def guard(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    module = load_script()
    monkeypatch.setattr(module.platform, "machine", lambda: "x86_64")
    return module


def test_sse4_listing_is_portable(guard: ModuleType) -> None:
    assert guard.find_non_portable_instructions(SSE4_LISTING) == []


@pytest.mark.parametrize(
    "listing, offending",
    [
        (AVX2_LISTING, "vpmaxsd %ymm1,%ymm2,%ymm0"),
        (AVX512_LISTING, "vpaddd %zmm1,%zmm2,%zmm0{%k1}"),
        (MASK_ONLY_LISTING, "kmovw  %k1,%eax"),
        (VEX_XMM_LISTING, "vmovdqa %xmm0,(%rdi)"),
    ],
)
def test_avx_family_instructions_are_flagged(
    guard: ModuleType, listing: str, offending: str
) -> None:
    assert guard.find_non_portable_instructions(listing) == [offending]


# -march=native also emits non-VEX-prefixed extensions beyond the SSE4.1
# baseline (POPCNT, LZCNT, BMI1/BMI2 general-purpose instructions, SSE4.2).
@pytest.mark.parametrize(
    "offending",
    [
        "popcnt %rdi,%rax",
        "lzcnt  %edi,%eax",
        "tzcnt  %rdi,%rax",
        "andn   %esi,%edi,%eax",
        "shlx   %rsi,%rdi,%rax",
        "sarx   %esi,%edi,%eax",
        "pext   %rsi,%rdi,%rax",
        "bzhi   %rsi,%rdi,%rax",
        "crc32q %rsi,%rax",
        "pcmpistri $0x0,%xmm1,%xmm0",
        "pcmpgtq %xmm1,%xmm0",
    ],
)
def test_non_baseline_scalar_extensions_are_flagged(guard: ModuleType, offending: str) -> None:
    listing = SSE4_LISTING + f"    1050:\t{offending}\n"
    assert guard.find_non_portable_instructions(listing) == [offending]


def test_sse41_lookalike_mnemonics_are_not_flagged(guard: ModuleType) -> None:
    # pextrd/pextrq are SSE4.1 (not BMI2 pext); pcmpeqq is SSE4.1 (not SSE4.2 pcmpgtq).
    extra = "".join(
        f"    {1060 + i:x}:\t{text}\n"
        for i, text in enumerate(
            ["pextrd $0x3,%xmm2,%eax", "pextrq $0x1,%xmm0,%rax", "pcmpeqq %xmm1,%xmm0"]
        )
    )
    assert guard.find_non_portable_instructions(SSE4_LISTING + extra) == []


def _fake_extension(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, guard: ModuleType) -> Path:
    extension = tmp_path / "pyabpoa.cpython-312-x86_64-linux-gnu.so"
    extension.write_bytes(b"\x7fELF")
    monkeypatch.setattr(guard, "locate_extension", lambda name: extension)
    return extension


def _fake_objdump(monkeypatch: pytest.MonkeyPatch, guard: ModuleType, listing: str) -> list[Any]:
    calls: list[Any] = []

    def run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, stdout=listing, stderr="")

    monkeypatch.setattr(guard.subprocess, "run", run)
    return calls


def test_main_passes_portable_extension(
    guard: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    extension = _fake_extension(tmp_path, monkeypatch, guard)
    calls = _fake_objdump(monkeypatch, guard, SSE4_LISTING)
    assert guard.main(["pyabpoa"]) == 0
    assert calls == [["objdump", "-d", "--no-show-raw-insn", str(extension)]]
    assert "pyabpoa: portable" in capsys.readouterr().out


def test_main_fails_native_extension_and_names_remedy(
    guard: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _fake_extension(tmp_path, monkeypatch, guard)
    _fake_objdump(monkeypatch, guard, AVX512_LISTING)
    assert guard.main(["pyabpoa", "--objdump", "/usr/bin/objdump"]) == 1
    out = capsys.readouterr().out
    assert "pyabpoa: 1 non-baseline instruction" in out
    assert "vpaddd %zmm1,%zmm2,%zmm0{%k1}" in out
    assert "uv cache clean pyabpoa" in out


def test_main_reports_missing_objdump(
    guard: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _fake_extension(tmp_path, monkeypatch, guard)

    def missing(command: list[str], **kwargs: Any) -> None:
        raise FileNotFoundError(command[0])

    monkeypatch.setattr(guard.subprocess, "run", missing)
    assert guard.main(["pyabpoa"]) == 2
    assert "objdump not found" in capsys.readouterr().err


def test_main_reports_failed_objdump(
    guard: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _fake_extension(tmp_path, monkeypatch, guard)

    def failing(command: list[str], **kwargs: Any) -> None:
        raise subprocess.CalledProcessError(1, command, stderr="bad file")

    monkeypatch.setattr(guard.subprocess, "run", failing)
    assert guard.main(["pyabpoa"]) == 2
    assert "bad file" in capsys.readouterr().err


def test_locate_extension_rejects_pure_python_module(guard: ModuleType) -> None:
    with pytest.raises(guard.GuardError, match="not a compiled extension"):
        guard.locate_extension("json")


def test_locate_extension_rejects_missing_module(guard: ModuleType) -> None:
    with pytest.raises(guard.GuardError, match="not installed"):
        guard.locate_extension("no_such_module_for_portable_guard")


def test_locate_extension_finds_installed_pyabpoa(guard: ModuleType) -> None:
    import pyabpoa

    assert guard.locate_extension("pyabpoa") == Path(pyabpoa.__file__)


def test_main_reports_guard_error(guard: ModuleType, capsys: pytest.CaptureFixture[str]) -> None:
    assert guard.main(["no_such_module_for_portable_guard"]) == 2
    assert "not installed" in capsys.readouterr().err


def test_main_skips_non_x86_64_hosts(
    guard: ModuleType, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(guard.platform, "machine", lambda: "aarch64")
    monkeypatch.setattr(guard, "locate_extension", pytest.fail)
    assert guard.main(["pyabpoa"]) == 0
    assert "skipped" in capsys.readouterr().out
