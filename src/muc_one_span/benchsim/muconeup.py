"""Argument builders for the external MucOneUp CLI (never imported as a library)."""

from __future__ import annotations

import re
from pathlib import Path

from muc_one_span.tools import run_tool

MIN_MUCONEUP = (0, 45, 0)
BUILTIN_PROFILE = {
    "ont_amplicon_r10": "ont_r10_sup_amplicon_v1",
    "ont_genomic_targeted": "ont_r10_genomic_v1",
    "hifi_amplicon": "hifi_amplicon_v1",
}


def require_muconeup(executable: str) -> str:
    """Return the MucOneUp version; fail if it predates read profiles and read truth."""
    out = run_tool([executable, "--version"])
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    if not match or tuple(map(int, match.groups())) < MIN_MUCONEUP:
        raise RuntimeError(f"MucOneUp >= 0.45.0 is required, found: {out.strip()!r}")
    return match.group(0)


def simulate_args(
    executable: str,
    config: Path,
    out_dir: Path,
    base: str,
    seed: int,
    lengths: tuple[int, int] | None,
    structure_file: Path | None,
    mutation: str | None,
    targets: tuple[tuple[int, int], ...],
) -> list[str]:
    """`muconeup simulate` for one diploid design."""
    args = [
        executable,
        "--config",
        str(config),
        "simulate",
        "--out-dir",
        str(out_dir),
        "--out-base",
        base,
        "--num-haplotypes",
        "2",
        "--output-structure",
        "--seed",
        str(seed),
    ]
    if structure_file is not None:
        args += ["--input-structure", str(structure_file)]
    elif lengths is not None:
        for length in lengths:
            args += ["--fixed-lengths", str(length)]
    else:
        raise ValueError("either lengths or structure_file is required")
    if mutation:
        args += ["--mutation-name", mutation]
        for hap, repeat in targets:
            args += ["--mutation-targets", f"{hap},{repeat}"]
    return args


def reads_args(
    executable: str,
    config: Path,
    profile: str,
    truth_fa: Path,
    out_dir: Path,
    base: str,
    seed: int,
    *,
    amount: int,
    profile_path: Path | None,
    flank_fasta: Path | None,
    pcr_preset: str | None,
) -> list[str]:
    """`muconeup reads …` for one benchmark profile; always FASTQ (`--no-align`)."""
    ref = str(profile_path) if profile_path else BUILTIN_PROFILE[profile]
    head = [executable, str(config), "--out-dir", str(out_dir)]
    tail = [
        "--read-profile",
        ref,
        "--out-base",
        base,
        "--seed",
        str(seed),
        "--no-align",
        str(truth_fa),
    ]
    if profile == "ont_genomic_targeted":
        args = [*head, "reads", "ont", "--simulator", "pbsim3-fragments", "--n-reads", str(amount)]
        if flank_fasta is not None:
            args.extend(["--flank-fasta", str(flank_fasta)])
        return [*args, *tail]
    platform = "ont" if profile == "ont_amplicon_r10" else "pacbio"
    args = [*head, "reads", "amplicon", "--platform", platform, "--coverage", str(amount)]
    if pcr_preset:
        args.extend(["--pcr-preset", pcr_preset])
    return [*args, *tail]
