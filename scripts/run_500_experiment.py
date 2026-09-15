#!/usr/bin/env python3
"""Execute 500-dataset MucOneUp experiment with single-truth generation and dual ledgers.

Supports Amplicon HiFi, Genomic ONT (NanoSim), and Genomic PacBio (PBSIM3/CCS).
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from muc_one_span.config import _bundled_repeats_path, load_repeat_dictionary  # noqa: E402
from muc_one_span.durable_ledger import DurableLedger, LedgerEntry, compute_sha256  # noqa: E402
from muc_one_span.evaluation.truth import load_truth  # noqa: E402
from muc_one_span.experiments import count_usable_records  # noqa: E402
from muc_one_span.inventory import (  # noqa: E402
    SampleInventoryRecord,
    build_expected_inventory,
)
from muc_one_span.tools import run_tool  # noqa: E402

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
        default=shutil.which("muconeup") or "/home/bernt-popp/miniforge3/bin/muconeup",
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
        choices=["primary", "amplicon_hifi", "genomic_ont", "genomic_pacbio", "all"],
        default="primary",
        help="Platforms to simulate ('primary' = amplicon_hifi + genomic_ont; 'all' adds genomic_pacbio)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of concurrent simulation workers (default: 4)",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional maximum number of samples to process",
    )
    parser.add_argument(
        "--unseal-test",
        action="store_true",
        help="Explicit confirmation to generate or access held-out test split",
    )
    parser.add_argument(
        "--categories",
        nargs="*",
        default=None,
        help="Optional list of categories to filter",
    )
    return parser


def simulate_biological_truth(
    record: SampleInventoryRecord,
    config_path: Path,
    truth_dir: Path,
    executable: str,
) -> Path:
    """Generate biological truth once per design using muconeup simulate."""
    truth_dir.mkdir(parents=True, exist_ok=True)
    out_base = record.design_name
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
        str(record.bio_seed),
    ]
    if record.structure_file:
        struct_path = ROOT / record.structure_file
        cmd.extend(["--input-structure", str(struct_path)])
    else:
        for length in record.lengths:
            cmd.extend(["--fixed-lengths", str(length)])

    if record.mutation and record.targets:
        mut_name = record.mutation
        if mut_name == "60dupA":
            mut_name = "dupA"
        cmd.extend(["--mutation-name", mut_name])
        for hap, rep in record.targets:
            cmd.extend(["--mutation-targets", f"{hap},{rep}"])

    run_tool(cmd, cwd=str(config_path.parent))
    if not expected_fa.exists() or expected_fa.stat().st_size == 0:
        raise FileNotFoundError(f"Expected simulated FASTA not created: {expected_fa}")
    return expected_fa


def _sanitize_fastq_headers(fastq_path: Path) -> None:
    """Ensure every read header has a unique deterministic record index suffix."""
    lines = fastq_path.read_text().splitlines()
    if not lines or len(lines) < 4 or "_r" in lines[0].split()[0]:
        return
    sanitized: list[str] = []
    rec_idx = 0
    for i in range(0, len(lines), 4):
        if i + 3 >= len(lines):
            break
        header = lines[i]
        seq = lines[i + 1]
        plus = lines[i + 2]
        qual = lines[i + 3]
        if header.startswith("@"):
            parts = header.split(None, 1)
            new_id = f"{parts[0]}_r{rec_idx:04d}"
            header = f"{new_id} {parts[1]}" if len(parts) > 1 else new_id
            rec_idx += 1
        sanitized.extend([header, seq, plus, qual])
    tmp_san = fastq_path.with_suffix(".san.tmp")
    tmp_san.write_text("\n".join(sanitized) + "\n", encoding="utf-8")
    with tmp_san.open("rb") as f:
        os.fsync(f.fileno())
    tmp_san.replace(fastq_path)


def simulate_platform_reads(
    record: SampleInventoryRecord,
    truth_fa: Path,
    config_path: Path,
    output_dir: Path,
    executable: str,
) -> Path:
    """Simulate platform reads for a given acquisition mode."""
    output_dir.mkdir(parents=True, exist_ok=True)
    platform = record.platform_mode
    out_base = f"{record.design_name}_{platform}"
    seed = record.platform_seed
    coverage = record.requested_templates

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
    if not candidates:
        raise FileNotFoundError(f"No read outputs generated in {output_dir}")
    read_file = candidates[0]
    if read_file.name.endswith((".fastq", ".fq")):
        _sanitize_fastq_headers(read_file)
    return read_file


def _worker_simulate(
    record: SampleInventoryRecord,
    truth_fa: Path,
    config_path: Path,
    output_dir: Path,
    executable: str,
) -> LedgerEntry:
    """Worker task to simulate platform reads and build structured LedgerEntry."""
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    plat_reads_dir = output_dir / "reads" / record.design_name / record.platform_mode
    reads_file = simulate_platform_reads(record, truth_fa, config_path, plat_reads_dir, executable)
    reads_hash = compute_sha256(reads_file)
    truth_hash = compute_sha256(truth_fa)
    usable_count = count_usable_records(reads_file)

    # Link verified atomic blinded copy
    blinded_dir = output_dir / "blinded_reads" / record.token / record.platform_mode
    blinded_dir.mkdir(parents=True, exist_ok=True)
    blinded_file = blinded_dir / f"{record.token}_{record.platform_mode}{reads_file.suffix}"
    if not blinded_file.exists() or compute_sha256(blinded_file) != reads_hash:
        tmp_blinded = blinded_file.with_suffix(".tmp")
        shutil.copy2(reads_file, tmp_blinded)
        with tmp_blinded.open("rb") as f:
            os.fsync(f.fileno())
        tmp_blinded.replace(blinded_file)
        if compute_sha256(blinded_file) != reads_hash:
            raise RuntimeError(f"Blinded copy hash mismatch for {blinded_file}")

    git_sha = os.environ.get("GIT_COMMIT_SHA")
    config_hash = compute_sha256(config_path)

    return LedgerEntry(
        design_id=record.design_id,
        design_name=record.design_name,
        token=record.token,
        split=record.split,
        category=record.category,
        platform=record.platform_mode,
        lengths=record.lengths,
        mutation=record.mutation,
        targets=record.targets,
        truth_fa=str(truth_fa.relative_to(output_dir)),
        truth_sha256=truth_hash,
        reads_file=str(reads_file.relative_to(output_dir)),
        reads_sha256=reads_hash,
        usable_records=usable_count,
        git_sha=git_sha,
        config_sha256=config_hash,
        status="completed",
    )


def run_pilot(config_path: Path, executable: str, output_dir: Path) -> None:
    """Run pilot validation across Amplicon, Genomic ONT, and Genomic PacBio."""
    pilot_dir = output_dir / "pilot_500"
    pilot_dir.mkdir(parents=True, exist_ok=True)
    print("=== Launching Pilot Run for 500-Dataset Pipeline ===")

    recs = [
        SampleInventoryRecord(
            sample_id="pilot_ctrl_60_60_amplicon_hifi",
            token="sample_9991",
            design_id=9991,
            design_name="pilot_ctrl_60_60",
            split="pilot",
            category="C1_HOMOZYGOUS_WT",
            platform_mode="amplicon_hifi",
            sequencing_platform="pacbio",
            lengths=[60, 60],
            mutation=None,
            targets=[],
            structure_file="examples/structures_500/identical_60_60.txt",
            bio_seed=999100,
            platform_seed=999110,
            requested_templates=100,
        ),
        SampleInventoryRecord(
            sample_id="pilot_mut_dupc_60_80_amplicon_hifi",
            token="sample_9992",
            design_id=9992,
            design_name="pilot_mut_dupc_60_80",
            split="pilot",
            category="C4_PATHOGENIC_DUPC",
            platform_mode="amplicon_hifi",
            sequencing_platform="pacbio",
            lengths=[60, 80],
            mutation="dupC",
            targets=[[2, 45]],
            structure_file=None,
            bio_seed=999200,
            platform_seed=999210,
            requested_templates=100,
        ),
    ]

    rd = load_repeat_dictionary(_bundled_repeats_path())
    for r in recs:
        print(f"Generating biological truth for {r.design_name}...")
        t_dir = pilot_dir / r.design_name / "truth"
        truth_fa = simulate_biological_truth(r, config_path, t_dir, executable)
        truth_hash = compute_sha256(truth_fa)
        truth_obj = load_truth(t_dir, rd)
        print(f"  Truth: {truth_hash[:16]}... Haplotypes: {len(truth_obj.haplotypes)}")

        for plat in ["amplicon_hifi", "genomic_ont", "genomic_pacbio"]:
            p_rec = SampleInventoryRecord(
                sample_id=f"{r.design_name}_{plat}",
                token=r.token,
                design_id=r.design_id,
                design_name=r.design_name,
                split="pilot",
                category=r.category,
                platform_mode=plat,
                sequencing_platform="ont" if plat == "genomic_ont" else "pacbio",
                lengths=r.lengths,
                mutation=r.mutation,
                targets=r.targets,
                structure_file=r.structure_file,
                bio_seed=r.bio_seed,
                platform_seed=r.platform_seed,
                requested_templates=r.requested_templates,
            )
            print(f"  Simulating {plat} reads...")
            reads_file = simulate_platform_reads(
                p_rec, truth_fa, config_path, pilot_dir / r.design_name / plat, executable
            )
            usable = count_usable_records(reads_file)
            print(f"  {plat} complete. Records: {usable}, Path: {reads_file.name}")

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

    if args.split in ("test", "all") and not args.unseal_test:
        print(
            "Error: Generating or accessing the held-out test split requires explicit --unseal-test",
            file=sys.stderr,
        )
        return 1

    include_supp = args.platforms in ("all", "genomic_pacbio")
    all_records = build_expected_inventory(designs_path, include_supplementary=include_supp)

    # Filter split
    if args.split != "all":
        selected = [r for r in all_records if r.split == args.split]
    else:
        selected = all_records

    # Filter platforms
    if args.platforms == "primary":
        selected = [r for r in selected if r.platform_mode in ("amplicon_hifi", "genomic_ont")]
    elif args.platforms != "all":
        selected = [r for r in selected if r.platform_mode == args.platforms]

    # Filter categories
    if args.categories:
        selected = [r for r in selected if r.category in args.categories]

    if args.max_samples:
        selected = selected[: args.max_samples]

    ledger = DurableLedger(output_dir)
    with ledger.run_lock():
        truth_root = output_dir / "truth"
        truth_root.mkdir(parents=True, exist_ok=True)

        print(
            f"=== Starting Generation for split='{args.split}', platforms='{args.platforms}' "
            f"({len(selected)} total datasets, workers={args.workers}) ==="
        )

        # Phase 1: Ensure biological truth exists for each distinct design
        unique_designs: dict[str, SampleInventoryRecord] = {}
        for r in selected:
            if r.design_name not in unique_designs:
                unique_designs[r.design_name] = r

        truth_paths: dict[str, Path] = {}
        for d_name, r in unique_designs.items():
            d_truth_dir = truth_root / d_name
            t_path = simulate_biological_truth(r, config_path, d_truth_dir, args.muconeup)
            truth_paths[d_name] = t_path

        # Phase 2: Parallel platform read simulation with durable atomic commits
        to_process: list[SampleInventoryRecord] = []
        skipped_count = 0
        for r in selected:
            if ledger.is_verified_complete(r.design_name, r.platform_mode):
                skipped_count += 1
            else:
                to_process.append(r)

        print(f"Already verified complete: {skipped_count}. To simulate: {len(to_process)}")

        if to_process:
            executor = ThreadPoolExecutor(max_workers=max(1, min(args.workers, 8)))
            future_to_rec = {
                executor.submit(
                    _worker_simulate,
                    r,
                    truth_paths[r.design_name],
                    config_path,
                    output_dir,
                    args.muconeup,
                ): r
                for r in to_process
            }
            completed = 0
            try:
                for fut in as_completed(future_to_rec):
                    rec = future_to_rec[fut]
                    entry = fut.result()
                    ledger.commit_entry(entry)
                    completed += 1
                    print(
                        f"[{completed}/{len(to_process)}] Completed and committed: "
                        f"{rec.design_name} ({rec.platform_mode})"
                    )
            except (KeyboardInterrupt, BaseException) as exc:
                print(f"Aborting generation batch: {exc}", file=sys.stderr)
                for f in future_to_rec:
                    f.cancel()
                executor.shutdown(wait=False, cancel_futures=True)
                raise
            finally:
                executor.shutdown(wait=True)

    print(f"Generation complete. Authoritative ledger at: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
