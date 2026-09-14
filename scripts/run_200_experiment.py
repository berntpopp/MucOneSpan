#!/usr/bin/env python3
"""Execute 200-dataset MucOneUp experiment with single-truth generation and dual ledgers."""

from __future__ import annotations

import argparse
import json
import random
import shutil
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from muc_one_span.config import _bundled_repeats_path, load_repeat_dictionary  # noqa: E402
from muc_one_span.evaluation.truth import load_truth  # noqa: E402
from muc_one_span.experiments import _hash, count_usable_records  # noqa: E402
from muc_one_span.tools import run_tool  # noqa: E402

# Locating generation_config.json
candidate_configs = [
    Path(
        "/home/bernt-popp/development/MucOneSpan/.worktrees/production-validation/tests/results/production_validation_20260914/fresh_generation/prepared_v2/generation_config.json"
    ),
    ROOT.parent
    / "production-validation"
    / "tests"
    / "results"
    / "production_validation_20260914"
    / "fresh_generation"
    / "prepared_v2"
    / "generation_config.json",
    ROOT
    / "tests"
    / "results"
    / "production_validation_20260914"
    / "fresh_generation"
    / "prepared_v2"
    / "generation_config.json",
]
DEFAULT_CONFIG = next((p for p in candidate_configs if p.exists()), candidate_configs[0])


def parse_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--designs",
        type=Path,
        default=ROOT / "examples" / "experiment_100_designs.json",
        help="Path to 100-design JSON specification",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Path to simulator configuration JSON",
    )
    parser.add_argument(
        "--muconeup",
        default="/home/bernt-popp/miniforge3/bin/muconeup",
        help="Path to muconeup executable",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "tests" / "data" / "experiment_200",
        help="Base output directory for generated datasets",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="Run pilot validation only (4 datasets, excluded from 200-case set)",
    )
    parser.add_argument(
        "--split",
        choices=["dev", "final", "all"],
        default="dev",
        help="Which split to generate (default: dev)",
    )
    return parser


def simulate_biological_truth(
    design: dict[str, Any],
    config_path: Path,
    truth_dir: Path,
    executable: str,
) -> Path:
    """Generate biological truth once using muconeup simulate."""
    truth_dir.mkdir(parents=True, exist_ok=True)
    out_base = design["name"]
    expected_fa = truth_dir / f"{out_base}.001.simulated.fa"
    if expected_fa.exists():
        return expected_fa

    cmd = [
        executable,
        "--config",
        str(config_path),
        "simulate",
        "--out-base",
        out_base,
        "--out-dir",
        str(truth_dir),
        "--num-haplotypes",
        "2",
        "--output-structure",
        "--seed",
        str(design["bio_seed"]),
    ]
    if design.get("structure_file"):
        struct_path = ROOT / design["structure_file"]
        cmd.extend(["--input-structure", str(struct_path)])
    else:
        for length in design["lengths"]:
            cmd.extend(["--fixed-lengths", str(length)])

    if design.get("mutation") and design.get("targets"):
        cmd.extend(["--mutation-name", design["mutation"]])
        for hap, rep in design["targets"]:
            cmd.extend(["--mutation-targets", f"{hap},{rep}"])

    run_tool(cmd, cwd=str(config_path.parent))
    if not expected_fa.exists():
        raise FileNotFoundError(f"Expected simulated FASTA not created: {expected_fa}")
    return expected_fa


def _sanitize_fastq_headers(fastq_path: Path) -> None:
    """Ensure every read header has a unique deterministic record index suffix."""
    lines = fastq_path.read_text().splitlines()
    if not lines or len(lines) < 4:
        return
    if "_r" in lines[0].split()[0]:
        return
    sanitized: list[str] = []
    rec_idx = 0
    for i in range(0, len(lines), 4):
        header = lines[i]
        if header.startswith("@"):
            parts = header.split(None, 1)
            new_id = f"{parts[0]}_r{rec_idx:04d}"
            header = f"{new_id} {parts[1]}" if len(parts) > 1 else new_id
            rec_idx += 1
        sanitized.append(header)
        sanitized.extend(lines[i + 1 : min(i + 4, len(lines))])
    fastq_path.write_text("\n".join(sanitized) + "\n")


def simulate_platform_reads(
    design: dict[str, Any],
    truth_fa: Path,
    platform: str,
    config_path: Path,
    output_dir: Path,
    executable: str,
) -> Path:
    """Simulate platform amplicon reads from the shared biological truth FASTA."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_base = f"{design['name']}_{platform}"
    seed = design["hifi_seed"] if platform == "hifi" else design["ont_seed"]
    sim_platform = "pacbio" if platform == "hifi" else "ont"

    expected_fastq = output_dir / f"{out_base}_reads_amplicon_{sim_platform}.fastq"
    if expected_fastq.exists() and expected_fastq.stat().st_size > 0:
        _sanitize_fastq_headers(expected_fastq)
        return expected_fastq

    cmd = [
        executable,
        "--config",
        str(config_path),
        "reads",
        "amplicon",
        str(truth_fa),
        "--out-dir",
        str(output_dir),
        "--out-base",
        out_base,
        "--coverage",
        str(design["requested_templates"]),
        "--seed",
        str(seed),
        "--platform",
        sim_platform,
    ]
    run_tool(cmd, cwd=str(config_path.parent))

    # Identify the generated read file (FASTQ)
    expected_fastq = output_dir / f"{out_base}_reads_amplicon_{sim_platform}.fastq"
    if not expected_fastq.exists():
        # Check alternative naming
        candidates = (
            list(output_dir.glob("*.fastq"))
            + list(output_dir.glob("*.fq"))
            + list(output_dir.glob("*.fq.gz"))
        )
        if len(candidates) == 1:
            expected_fastq = candidates[0]
        else:
            raise FileNotFoundError(
                f"Could not locate unique simulated read output in {output_dir}; candidates={candidates}"
            )
    _sanitize_fastq_headers(expected_fastq)
    return expected_fastq


def run_pilot(config_path: Path, executable: str, output_dir: Path) -> None:
    """Run a 2-design x 2-platform pilot to validate generation mechanics."""
    pilot_dir = output_dir / "pilot"
    pilot_dir.mkdir(parents=True, exist_ok=True)
    print("=== Launching Pilot Run (4 datasets) ===")

    pilot_designs = [
        {
            "design_id": 9991,
            "name": "pilot_ctrl_60_60",
            "split": "pilot",
            "lengths": [60, 60],
            "mutation": None,
            "targets": [],
            "requested_templates": 200,
            "structure_file": "examples/structures/identical_60_60.txt",
            "bio_seed": 999100,
            "hifi_seed": 999110,
            "ont_seed": 999120,
        },
        {
            "design_id": 9992,
            "name": "pilot_mut_asym_25_140",
            "split": "pilot",
            "lengths": [25, 140],
            "mutation": "dupC",
            "targets": [[2, 100]],
            "requested_templates": 60,
            "structure_file": None,
            "bio_seed": 999200,
            "hifi_seed": 999210,
            "ont_seed": 999220,
        },
    ]

    rd = load_repeat_dictionary(_bundled_repeats_path())
    for d in pilot_designs:
        d_name = str(d["name"])
        print(f"Generating biological truth for {d_name}...")
        truth_fa = simulate_biological_truth(
            d, config_path, pilot_dir / d_name / "truth", executable
        )
        truth_hash = _hash(truth_fa)
        truth_obj = load_truth(pilot_dir / d_name / "truth", rd)
        print(
            f"  Truth generated successfully. SHA256: {truth_hash[:16]}... Haplotypes: {len(truth_obj.haplotypes)}"
        )

        for p in ["hifi", "ont"]:
            print(f"  Simulating {p.upper()} reads...")
            reads_path = simulate_platform_reads(
                d, truth_fa, p, config_path, pilot_dir / d_name / p, executable
            )
            usable = count_usable_records(reads_path)
            print(f"  {p.upper()} complete. Usable records: {usable}, Path: {reads_path.name}")

    print("=== Pilot Run Completed Successfully ===")


def main() -> int:
    args = parse_args().parse_args()
    if not args.config.exists():
        print(f"Error: Simulator config not found: {args.config}", file=sys.stderr)
        return 1

    if args.pilot:
        run_pilot(args.config, args.muconeup, args.output_dir)
        return 0

    if not args.designs.exists():
        print(f"Error: Designs file not found: {args.designs}", file=sys.stderr)
        return 1

    with args.designs.open() as f:
        design_data = json.load(f)
    all_designs = design_data["designs"]

    # Filter by requested split
    if args.split == "dev":
        selected_designs = [d for d in all_designs if d["split"] == "dev"]
    elif args.split == "final":
        selected_designs = [d for d in all_designs if d["split"] == "final"]
    else:
        selected_designs = all_designs

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    truth_root = output_dir / "truth"
    reads_root = output_dir / "reads"
    blinded_root = output_dir / "blinded_reads"
    truth_root.mkdir(exist_ok=True)
    reads_root.mkdir(exist_ok=True)
    blinded_root.mkdir(exist_ok=True)

    public_ledger_path = output_dir / "ledger_public.jsonl"
    sealed_ledger_path = output_dir / "ledger_sealed.jsonl"

    rd = load_repeat_dictionary(_bundled_repeats_path())
    rng = random.Random(42_2026_0914)  # Seeded token shuffle

    print(
        f"=== Starting Generation for split='{args.split}': {len(selected_designs)} designs x 2 platforms = {len(selected_designs) * 2} datasets ==="
    )

    # Assign randomized neutral tokens for blinded caller directory
    token_map = {}
    tokens = [f"sample_{i:04d}" for i in range(1, len(all_designs) + 1)]
    rng.shuffle(tokens)
    for d, tok in zip(all_designs, tokens, strict=True):
        token_map[d["name"]] = tok

    existing_sealed_keys: set[tuple[str, str]] = set()
    if sealed_ledger_path.exists():
        with sealed_ledger_path.open() as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        e = json.loads(line)
                        existing_sealed_keys.add((e["design_name"], e["platform"]))
                    except Exception:
                        pass

    for idx, d in enumerate(selected_designs, 1):
        print(
            f"[{idx}/{len(selected_designs)}] Processing design '{d['name']}' (Split: {d['split']})..."
        )
        d_truth_dir = truth_root / d["name"]
        truth_fa = simulate_biological_truth(d, args.config, d_truth_dir, args.muconeup)
        truth_hash = _hash(truth_fa)

        # Validate truth
        truth_obj = load_truth(d_truth_dir, rd)
        if len(truth_obj.haplotypes) != 2:
            raise ValueError(
                f"Design {d['name']} generated {len(truth_obj.haplotypes)} haplotypes, expected 2"
            )

        for platform in ["hifi", "ont"]:
            if (d["name"], platform) in existing_sealed_keys:
                token = token_map[d["name"]]
                print(
                    f"  {platform.upper()} already in ledger; skipping (Token: {token}_{platform})"
                )
                continue
            p_output_dir = reads_root / d["name"] / platform
            p_start = time.monotonic()
            reads_file = simulate_platform_reads(
                d, truth_fa, platform, args.config, p_output_dir, args.muconeup
            )
            p_time = time.monotonic() - p_start
            reads_hash = _hash(reads_file)
            usable_count = count_usable_records(reads_file)

            # Blinded link / copy
            token = token_map[d["name"]]
            blinded_name = f"{token}_{platform}.fastq"
            blinded_path = blinded_root / blinded_name
            if not blinded_path.exists():
                shutil.copyfile(reads_file, blinded_path)

            # Public ledger entry
            pub_entry = {
                "sample_token": token,
                "platform": platform,
                "split": d["split"],
                "input_fastq": str(blinded_path.relative_to(ROOT)),
                "input_sha256": reads_hash,
                "usable_records": usable_count,
                "status": "completed",
                "wall_seconds": round(p_time, 2),
            }
            with public_ledger_path.open("a") as f_pub:
                f_pub.write(json.dumps(pub_entry) + "\n")

            # Sealed truth ledger entry
            sealed_entry = {
                "sample_token": token,
                "design_id": d["design_id"],
                "design_name": d["name"],
                "split": d["split"],
                "platform": platform,
                "category": d["category"],
                "lengths": d["lengths"],
                "mutation": d["mutation"],
                "targets": d["targets"],
                "bio_seed": d["bio_seed"],
                "platform_seed": d["hifi_seed"] if platform == "hifi" else d["ont_seed"],
                "truth_fa_sha256": truth_hash,
                "raw_reads_path": str(reads_file.relative_to(ROOT)),
                "raw_reads_sha256": reads_hash,
                "usable_records": usable_count,
            }
            with sealed_ledger_path.open("a") as f_seal:
                f_seal.write(json.dumps(sealed_entry) + "\n")

            print(
                f"  {platform.upper()} done: {usable_count} usable reads in {p_time:.1f}s (Token: {token}_{platform})"
            )

    print(f"\nGeneration complete for split='{args.split}'.")
    print(f"Public ledger: {public_ledger_path}")
    print(f"Sealed ledger: {sealed_ledger_path} (SHA-256: {_hash(sealed_ledger_path)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
