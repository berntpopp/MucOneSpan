#!/usr/bin/env python3
"""Execute 500-dataset MucOneUp experiment with single-truth generation and dual ledgers.

Supports Amplicon HiFi, Genomic ONT (Adaptive Sampling / NanoSim), and Genomic PacBio (PBSIM3/CCS).
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from muc_one_span.config import _bundled_repeats_path, load_repeat_dictionary  # noqa: E402
from muc_one_span.evaluation.truth import load_truth  # noqa: E402
from muc_one_span.experiments import _hash, count_usable_records  # noqa: E402
from muc_one_span.tools import run_tool  # noqa: E402

# Candidate paths for generation_config.json
candidate_configs = [
    ROOT / "examples" / "generation_config.json",
    ROOT.parent / "MucOneUp" / "config" / "generation_config.json",
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
        default=ROOT / "examples" / "experiment_500_designs.json",
        help="Path to 250-design JSON specification",
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
        default=ROOT / "tests" / "data" / "experiment_500",
        help="Base output directory for generated datasets",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="Run pilot validation only (4 datasets across all platforms)",
    )
    parser.add_argument(
        "--split",
        choices=["dev", "val", "test", "all"],
        default="dev",
        help="Which split to generate (default: dev)",
    )
    parser.add_argument(
        "--platforms",
        choices=["amplicon_hifi", "genomic_ont", "genomic_pacbio", "all"],
        default="amplicon_hifi",
        help="Which sequencing mode / platform to simulate (default: amplicon_hifi)",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional maximum number of samples to process",
    )
    parser.add_argument(
        "--samples-per-category",
        type=int,
        default=None,
        help="Optional maximum number of samples to select per category",
    )
    parser.add_argument(
        "--categories",
        nargs="*",
        default=None,
        help="Optional list of categories to filter",
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
    if expected_fa.exists() and expected_fa.stat().st_size > 0:
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
        mut_name = design["mutation"]
        if mut_name == "60dupA":
            mut_name = "dupA"
        cmd.extend(["--mutation-name", mut_name])
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
    """Simulate platform reads (amplicon, genomic ONT, or genomic PacBio)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_base = f"{design['name']}_{platform}"
    seed = design["hifi_seed"] if "hifi" in platform or "pacbio" in platform else design["ont_seed"]
    coverage = design["requested_templates"]

    candidates = [
        p
        for p in (
            list(output_dir.glob("*.fastq"))
            + list(output_dir.glob("*.fq"))
            + list(output_dir.glob("*.fq.gz"))
            + list(output_dir.glob("*.bam"))
        )
        if p.is_file() and p.stat().st_size > 0
    ]
    if candidates:
        read_file = candidates[0]
        if read_file.name.endswith((".fastq", ".fq")):
            _sanitize_fastq_headers(read_file)
        return read_file

    if platform == "amplicon_hifi":
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
            str(coverage),
            "--seed",
            str(seed),
            "--platform",
            "pacbio",
        ]
    elif platform == "genomic_ont":
        cmd = [
            executable,
            "--config",
            str(config_path),
            "reads",
            "ont",
            str(truth_fa),
            "--out-dir",
            str(output_dir),
            "--out-base",
            out_base,
            "--coverage",
            str(max(10, coverage // 5)),
            "--seed",
            str(seed),
        ]
    elif platform == "genomic_pacbio":
        cmd = [
            executable,
            "--config",
            str(config_path),
            "reads",
            "pacbio",
            str(truth_fa),
            "--out-dir",
            str(output_dir),
            "--out-base",
            out_base,
            "--coverage",
            str(max(10, coverage // 5)),
            "--seed",
            str(seed),
        ]
    else:
        raise ValueError(f"Unsupported platform mode: {platform}")

    run_tool(cmd, cwd=str(config_path.parent))

    candidates = (
        list(output_dir.glob("*.fastq"))
        + list(output_dir.glob("*.fq"))
        + list(output_dir.glob("*.fq.gz"))
        + list(output_dir.glob("*.bam"))
    )
    if not candidates:
        raise FileNotFoundError(f"No read outputs generated in {output_dir}")
    read_file = candidates[0]
    if read_file.name.endswith((".fastq", ".fq")):
        _sanitize_fastq_headers(read_file)
    return read_file


def run_pilot(config_path: Path, executable: str, output_dir: Path) -> None:
    """Run pilot validation across Amplicon, Genomic ONT, and Genomic PacBio."""
    pilot_dir = output_dir / "pilot_500"
    pilot_dir.mkdir(parents=True, exist_ok=True)
    print("=== Launching Pilot Run for 500-Dataset Pipeline ===")

    pilot_designs = [
        {
            "design_id": 9991,
            "name": "pilot_ctrl_60_60",
            "split": "pilot",
            "lengths": [60, 60],
            "mutation": None,
            "targets": [],
            "requested_templates": 100,
            "structure_file": "examples/structures_500/identical_60_60.txt",
            "bio_seed": 999100,
            "hifi_seed": 999110,
            "ont_seed": 999120,
        },
        {
            "design_id": 9992,
            "name": "pilot_mut_dupc_60_80",
            "split": "pilot",
            "lengths": [60, 80],
            "mutation": "dupC",
            "targets": [[2, 45]],
            "requested_templates": 100,
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
            f"  Truth generated. SHA256: {truth_hash[:16]}... Haplotypes: {len(truth_obj.haplotypes)}"
        )

        for plat in ["amplicon_hifi", "genomic_ont", "genomic_pacbio"]:
            print(f"  Simulating {plat} reads...")
            reads_path = simulate_platform_reads(
                d, truth_fa, plat, config_path, pilot_dir / d_name / plat, executable
            )
            usable = count_usable_records(reads_path)
            print(f"  {plat} complete. Records: {usable}, Path: {reads_path.name}")

    print("=== Pilot Run Completed Successfully ===")


def main() -> int:
    args = parse_args().parse_args()
    config_path = args.config.resolve()
    output_dir = args.output_dir.resolve()
    designs_path = args.designs.resolve()

    if not config_path.exists():
        print(f"Error: Simulator config not found: {config_path}", file=sys.stderr)
        return 1

    if args.pilot:
        run_pilot(config_path, args.muconeup, output_dir)
        return 0

    if not designs_path.exists():
        print(f"Error: Designs file not found: {designs_path}", file=sys.stderr)
        return 1

    with designs_path.open() as f:
        all_designs = json.load(f)

    if args.split != "all":
        selected_designs = [d for d in all_designs if d["split"] == args.split]
    else:
        selected_designs = all_designs

    if args.categories:
        selected_designs = [d for d in selected_designs if d.get("category") in args.categories]

    if args.samples_per_category:
        by_cat: dict[str, list[dict[str, Any]]] = {}
        for d in selected_designs:
            by_cat.setdefault(d.get("category", "unknown"), []).append(d)
        selected_designs = []
        for cat in sorted(by_cat.keys()):
            cat_list = by_cat[cat]
            if cat == "C2_ASYMMETRIC_LENGTH":
                special = [d for d in cat_list if d.get("name") == "c2_asym_25_140_dev"]
                others = [d for d in cat_list if d.get("name") != "c2_asym_25_140_dev"]
                combined = special + others[: max(0, args.samples_per_category - len(special))]
                selected_designs.extend(combined)
            else:
                selected_designs.extend(cat_list[: args.samples_per_category])

    if args.max_samples:
        selected_designs = selected_designs[: args.max_samples]

    platforms = (
        ["amplicon_hifi", "genomic_ont", "genomic_pacbio"]
        if args.platforms == "all"
        else [args.platforms]
    )

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

    existing_sealed_keys: set[tuple[str, str]] = set()
    if sealed_ledger_path.exists():
        for line in sealed_ledger_path.read_text().splitlines():
            if line.strip():
                try:
                    ent = json.loads(line)
                    existing_sealed_keys.add(
                        (str(ent.get("design_name")), str(ent.get("platform")))
                    )
                except Exception:
                    pass

    rd = load_repeat_dictionary(_bundled_repeats_path())
    rng = random.Random(500_2026_0915)

    token_map = {}
    tokens = [f"sample_{i:04d}" for i in range(1, len(all_designs) + 1)]
    rng.shuffle(tokens)
    for d, tok in zip(all_designs, tokens, strict=True):
        token_map[d["name"]] = tok

    print(
        f"=== Starting Generation for split='{args.split}': "
        f"{len(selected_designs)} designs x {len(platforms)} platforms = "
        f"{len(selected_designs) * len(platforms)} datasets ==="
    )

    for idx, d in enumerate(selected_designs, 1):
        print(
            f"[{idx}/{len(selected_designs)}] Processing design '{d['name']}' ({d['category']})..."
        )
        d_truth_dir = truth_root / d["name"]
        truth_fa = simulate_biological_truth(d, args.config, d_truth_dir, args.muconeup)
        truth_hash = _hash(truth_fa)
        truth_obj = load_truth(d_truth_dir, rd)
        if len(truth_obj.haplotypes) != 2:
            raise ValueError(f"Expected 2 truth haplotypes, found {len(truth_obj.haplotypes)}")

        token = token_map[d["name"]]

        for plat in platforms:
            plat_reads_dir = reads_root / d["name"] / plat
            reads_file = simulate_platform_reads(
                d, truth_fa, plat, args.config, plat_reads_dir, args.muconeup
            )
            reads_hash = _hash(reads_file)
            usable_count = count_usable_records(reads_file)

            # Link blinded read copy
            blinded_plat_dir = blinded_root / token / plat
            blinded_plat_dir.mkdir(parents=True, exist_ok=True)
            blinded_file = blinded_plat_dir / f"{token}_{plat}{reads_file.suffix}"
            if not blinded_file.exists():
                shutil.copy2(reads_file, blinded_file)

            sealed_entry = {
                "design_id": d["design_id"],
                "design_name": d["name"],
                "token": token,
                "split": d["split"],
                "category": d["category"],
                "platform": plat,
                "lengths": d["lengths"],
                "mutation": d.get("mutation"),
                "targets": d.get("targets", []),
                "truth_fa": str(truth_fa.relative_to(output_dir)),
                "truth_sha256": truth_hash,
                "reads_file": str(reads_file.relative_to(output_dir)),
                "reads_sha256": reads_hash,
                "usable_records": usable_count,
            }
            if (d["name"], plat) not in existing_sealed_keys:
                with sealed_ledger_path.open("a", encoding="utf-8") as f:
                    f.write(json.dumps(sealed_entry) + "\n")

                public_entry = {
                    "token": token,
                    "platform": plat,
                    "reads_file": str(blinded_file.relative_to(output_dir)),
                    "usable_records": usable_count,
                }
                with public_ledger_path.open("a", encoding="utf-8") as f:
                    f.write(json.dumps(public_entry) + "\n")
                existing_sealed_keys.add((d["name"], plat))

    print(f"Generation complete. Ledgers written to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
