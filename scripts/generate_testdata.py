#!/usr/bin/env python3
"""Generate test data using MucOneUp for muconespan integration tests.

Requires:
  - MucOneUp installed and on PATH (pip install muc_one_up)
  - PacBio conda env activated (conda activate env_pacbio)
  - MucOneUp config.json available

Usage:
  python scripts/generate_testdata.py
  # Or via Makefile:
  make generate-testdata
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

# Output directory for generated test data
OUTPUT_DIR = Path("tests/data/generated")

# Path to MucOneUp config.json (sibling repo, relative to project root)
MUCONEUP_CONFIG = Path(__file__).parent.parent.parent / "MucOneUp" / "config.json"

# Test sample definitions
# (name, hap1_length, hap2_length, mutation_name, mutation_targets, seed)
SAMPLES = [
    # --- Core mutation types (60/80 baseline) ---
    ("sample_dupc_60_80", 60, 80, "dupC", "1,25", 1001),
    ("sample_dupa_60_80", 60, 80, "dupA", "1,25", 1002),
    ("sample_insg_60_80", 60, 80, "insG", "1,25", 1003),
    ("sample_dupcccc_60_80", 60, 80, "insCCCC", "1,25", 1004),
    ("sample_del_60_80", 60, 80, "del18_31", "1,25", 1005),
    ("sample_normal_60_80", 60, 80, "normal", "", 1006),
    ("sample_homozygous_60_60", 60, 60, "dupC", "1,25", 1007),
    ("sample_asymmetric_25_140", 25, 140, "dupC", "1,10", 1008),
    ("sample_short_25_30", 25, 30, "dupC", "1,10", 1009),
    ("sample_long_120_140", 120, 140, "dupC", "1,50", 1010),
    # --- QUAL tuning: dupC at multiple seeds (QUAL variance) ---
    ("sample_dupc_60_80_s2", 60, 80, "dupC", "1,25", 2001),
    ("sample_dupc_60_80_s3", 60, 80, "dupC", "1,25", 3001),
    ("sample_dupc_60_80_s4", 60, 80, "dupC", "1,25", 4001),
    ("sample_dupc_60_80_s5", 60, 80, "dupC", "1,25", 5001),
    # --- QUAL tuning: dupC at different allele lengths ---
    ("sample_dupc_40_50", 40, 50, "dupC", "1,15", 1011),
    ("sample_dupc_80_100", 80, 100, "dupC", "1,35", 1012),
    ("sample_dupc_100_120", 100, 120, "dupC", "1,45", 1013),
    # --- Close allele pairs (indel-valley splitting stress) ---
    ("sample_dupc_50_55", 50, 55, "dupC", "1,20", 1014),
    ("sample_dupc_50_57", 50, 57, "dupC", "1,20", 1015),
    ("sample_dupc_50_60", 50, 60, "dupC", "1,20", 1016),
    ("sample_normal_50_55", 50, 55, "normal", "", 1017),
    ("sample_normal_50_60", 50, 60, "normal", "", 1018),
    # --- Different mutations at long allele lengths (boundary artifact risk) ---
    ("sample_dupa_100_120", 100, 120, "dupA", "1,45", 1019),
    ("sample_insg_100_120", 100, 120, "insG", "1,45", 1020),
    ("sample_normal_100_120", 100, 120, "normal", "", 1021),
    # --- Lower coverage (50x) for QUAL sensitivity ---
    ("sample_dupc_60_80_cov50", 60, 80, "dupC", "1,25", 1022),
]

# ONT samples use the same haplotypes but with ONT read simulation
ONT_SAMPLES = [
    ("sample_ont_dupc_60_80", 60, 80, "dupC", "1,25", 1001),
    ("sample_ont_normal_60_80", 60, 80, "normal", "", 1006),
    ("sample_ont_dupa_60_80", 60, 80, "dupA", "1,25", 1002),
]

COVERAGE = 200

# Samples that use non-default coverage
LOW_COVERAGE_SAMPLES = {"sample_dupc_60_80_cov50": 50}


def run(cmd: list[str], desc: str) -> None:
    """Run a command and exit on failure."""
    print(f"  {desc}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  FAILED: {' '.join(cmd)}")
        print(f"  stderr: {result.stderr}")
        sys.exit(1)


def generate_sample(
    name: str,
    hap1_length: int,
    hap2_length: int,
    mutation_name: str,
    mutation_targets: str,
    seed: int,
    config_path: Path,
    output_dir: Path,
    coverage: int = COVERAGE,
) -> None:
    """Generate a single test sample using MucOneUp."""
    sample_dir = output_dir / name
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Simulate haplotypes
    sim_cmd = [
        "muconeup",
        "--config",
        str(config_path),
        "simulate",
        "--out-base",
        name,
        "--out-dir",
        str(sample_dir),
        "--num-haplotypes",
        "2",
        # MucOneUp v0.44+ uses Click's multiple=True: pass --fixed-lengths once per value
        "--fixed-lengths",
        str(hap1_length),
        "--fixed-lengths",
        str(hap2_length),
        "--output-structure",
        "--seed",
        str(seed),
    ]

    # For "normal" mutation, do not pass --mutation-name or --mutation-targets
    if mutation_name != "normal":
        sim_cmd.extend(["--mutation-name", mutation_name])
        if mutation_targets:
            sim_cmd.extend(["--mutation-targets", mutation_targets])

    run(sim_cmd, f"Simulating haplotypes for {name}")

    # Step 2: Generate PacBio HiFi amplicon reads
    # Find the generated FASTA files
    fastas = sorted(sample_dir.glob(f"{name}.*.simulated.fa"))
    if not fastas:
        print(f"  ERROR: No FASTA files found for {name}")
        sys.exit(1)

    reads_cmd = [
        "muconeup",
        "--config",
        str(config_path),
        "reads",
        "amplicon",
        *[str(f) for f in fastas],
        "--out-dir",
        str(sample_dir),
        "--out-base",
        f"{name}_reads",
        "--coverage",
        str(coverage),
        "--seed",
        str(seed),
        "--platform",
        "pacbio",
    ]

    run(reads_cmd, f"Generating PacBio HiFi reads for {name} ({coverage}x)")


def generate_ont_sample(
    name: str,
    hap1_length: int,
    hap2_length: int,
    mutation_name: str,
    mutation_targets: str,
    seed: int,
    config_path: Path,
    output_dir: Path,
    coverage: int = COVERAGE,
) -> None:
    """Generate a single ONT test sample using MucOneUp."""
    sample_dir = output_dir / name
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Simulate haplotypes (same as PacBio)
    sim_cmd = [
        "muconeup",
        "--config",
        str(config_path),
        "simulate",
        "--out-base",
        name,
        "--out-dir",
        str(sample_dir),
        "--num-haplotypes",
        "2",
        "--fixed-lengths",
        str(hap1_length),
        "--fixed-lengths",
        str(hap2_length),
        "--output-structure",
        "--seed",
        str(seed),
    ]

    if mutation_name != "normal":
        sim_cmd.extend(["--mutation-name", mutation_name])
        if mutation_targets:
            sim_cmd.extend(["--mutation-targets", mutation_targets])

    run(sim_cmd, f"Simulating haplotypes for {name}")

    # Step 2: Generate ONT amplicon reads
    fastas = sorted(sample_dir.glob(f"{name}.*.simulated.fa"))
    if not fastas:
        print(f"  ERROR: No FASTA files found for {name}")
        sys.exit(1)

    reads_cmd = [
        "muconeup",
        "--config",
        str(config_path),
        "reads",
        "amplicon",
        *[str(f) for f in fastas],
        "--out-dir",
        str(sample_dir),
        "--out-base",
        f"{name}_reads",
        "--coverage",
        str(coverage),
        "--seed",
        str(seed),
        "--platform",
        "ont",
    ]

    run(reads_cmd, f"Generating ONT reads for {name} ({coverage}x)")


def main() -> None:
    """Generate all test samples."""
    import os

    # Find MucOneUp config -- prefer env var, then sibling-repo default, then home fallback
    config_env = os.environ.get("MUCONEUP_CONFIG")
    if config_env:
        config = Path(config_env)
    else:
        config = MUCONEUP_CONFIG
        if not config.exists():
            alt = Path.home() / "development" / "MucOneUp" / "config.json"
            if alt.exists():
                config = alt
            else:
                print(f"ERROR: MucOneUp config.json not found at {config}")
                print("Set MUCONEUP_CONFIG env var or place config.json in expected location.")
                sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Generating {len(SAMPLES)} test samples in {OUTPUT_DIR}")
    print(f"Using MucOneUp config: {config}")
    print(f"Coverage: {COVERAGE}x PacBio HiFi amplicon")
    print()

    for i, (name, h1, h2, mut, targets, seed) in enumerate(SAMPLES, start=1):
        cov = LOW_COVERAGE_SAMPLES.get(name, COVERAGE)
        print(f"[{i}/{len(SAMPLES)}] {name} ({cov}x)")
        generate_sample(name, h1, h2, mut, targets, seed, config, OUTPUT_DIR, coverage=cov)
        print()

    print(f"Done. {len(SAMPLES)} PacBio samples generated in {OUTPUT_DIR}")

    # ONT samples
    print()
    print(f"Generating {len(ONT_SAMPLES)} ONT test samples in {OUTPUT_DIR}")
    print(f"Coverage: {COVERAGE}x ONT amplicon")
    print()

    for i, (name, h1, h2, mut, targets, seed) in enumerate(ONT_SAMPLES, start=1):
        print(f"[{i}/{len(ONT_SAMPLES)}] {name} ({COVERAGE}x)")
        generate_ont_sample(name, h1, h2, mut, targets, seed, config, OUTPUT_DIR)
        print()

    print(f"Done. {len(ONT_SAMPLES)} ONT samples generated in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
