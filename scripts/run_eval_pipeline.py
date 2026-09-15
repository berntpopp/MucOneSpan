#!/usr/bin/env python3
"""Run and evaluate MucOneSpan across experiment datasets with failure atlas generation."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from muc_one_span.config import _bundled_repeats_path, load_repeat_dictionary  # noqa: E402
from muc_one_span.evaluation import (  # noqa: E402
    aggregate,
    evaluate_sample,
    load_observation,
    load_truth,
)

MODELS = {
    "hifi": "/home/bernt-popp/miniforge3/envs/env_clair3/bin/models/hifi",
    "ont": "/home/bernt-popp/miniforge3/envs/env_clair3/bin/models/ont",
}


def parse_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ledger",
        type=Path,
        default=ROOT / "tests" / "data" / "experiment_200" / "ledger_sealed.jsonl",
        help="Path to sealed ledger with truth mapping",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "tests" / "data" / "experiment_200",
        help="Base experiment data directory",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for pipeline and evaluation results",
    )
    parser.add_argument(
        "--split",
        choices=["dev", "final", "pilot"],
        default="dev",
        help="Which split to run (default: dev)",
    )
    parser.add_argument(
        "--platform",
        choices=["hifi", "ont", "amplicon_hifi", "genomic_ont", "genomic_pacbio", "all"],
        default="all",
        help="Filter by platform",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip pipeline execution if summary.json already exists",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Generate HTML report with embedded IGV for each sample",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads per sample run",
    )
    parser.add_argument(
        "--parallel-samples",
        type=int,
        default=1,
        help="Number of samples to run in parallel (default: 1)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional runtime config JSON for MucOneSpan",
    )
    parser.add_argument(
        "--src-dir",
        type=Path,
        default=None,
        help="Optional path to muc_one_span src directory (e.g. for baseline testing)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit to first N samples for testing",
    )
    return parser


def diagnose_first_failure(
    sample_res_dir: Path,
    truth_obj: Any,
    eval_row: dict[str, Any],
) -> dict[str, Any]:
    """Diagnose the FIRST incorrect pipeline stage for a failed or inaccurate sample."""
    # Stage 1: Usable reads
    status_file = sample_res_dir / "run_status.json"
    if not status_file.exists():
        return {"first_failure_stage": "execution", "reason": "run_status.json missing"}
    status_data = json.loads(status_file.read_text())
    status = status_data.get("status")
    if status == "insufficient_evidence":
        return {"first_failure_stage": "usable_reads", "reason": "insufficient_evidence"}
    if status == "execution_failed":
        return {"first_failure_stage": "execution", "reason": status_data.get("error")}

    # Stage 2: Initial Mapping
    mapping_bam = sample_res_dir / "mapping.bam"
    if not mapping_bam.exists():
        return {"first_failure_stage": "initial_mapping", "reason": "mapping.bam missing"}

    # Stage 3: Length Inference / Allele detection
    alleles_file = sample_res_dir / "alleles.json"
    if not alleles_file.exists():
        return {"first_failure_stage": "length_inference", "reason": "alleles.json missing"}
    alleles_data = json.loads(alleles_file.read_text())
    detected_lengths = sorted(
        [
            alleles_data[k]["length"]
            for k in ("allele_1", "allele_2")
            if k in alleles_data
            and isinstance(alleles_data[k], dict)
            and "length" in alleles_data[k]
        ]
    )
    true_lengths = sorted([len(h.structure) for h in truth_obj.haplotypes])

    # Check length correctness
    if eval_row.get("length_correct") is False or detected_lengths != true_lengths:
        return {
            "first_failure_stage": "length_inference",
            "reason": f"lengths mismatch: detected {detected_lengths} vs true {true_lengths}",
            "detected_lengths": detected_lengths,
            "true_lengths": true_lengths,
        }

    # Stage 5: Consensus & Trimming
    consensus_fastas = list(sample_res_dir.glob("consensus_allele_*.fa"))
    if not consensus_fastas:
        return {"first_failure_stage": "consensus", "reason": "consensus_allele_*.fa missing"}

    # Stage 6: Classification & Mutation Detection
    summary_file = sample_res_dir / "summary.json"
    if not summary_file.exists():
        return {"first_failure_stage": "classification", "reason": "summary.json missing"}

    metrics = eval_row.get("metrics", {})
    all_sequences_exact = metrics.get("all_sequences_exact", {}).get("min", 0) == 1
    event_tp = metrics.get("event_tp", {}).get("min", 0)
    event_fp = metrics.get("event_fp", {}).get("max", 0)
    total_truth_events = len([e for h in truth_obj.haplotypes for e in h.events])

    if not all_sequences_exact:
        seq_exact = metrics.get("sequence_exact", {}).get("min", 0)
        status = eval_row.get("status")
        if status == "ambiguous_reconstruction" and seq_exact >= 1 and event_fp == 0:
            return {"first_failure_stage": "none", "reason": "homozygous unphased sequence exact"}
        return {
            "first_failure_stage": "consensus",
            "reason": "exact sequence discordance despite correct length pair",
        }

    if event_tp < total_truth_events:
        return {
            "first_failure_stage": "variant_calling",
            "reason": f"mutation missed: tp={event_tp} vs true={total_truth_events}",
        }

    if event_fp > 0:
        return {
            "first_failure_stage": "variant_calling",
            "reason": f"extra false positive mutation call: fp={event_fp}",
        }

    return {"first_failure_stage": "none", "reason": "sample reconstructed correctly"}


def run_sample(
    sample_id: str,
    fastq_path: Path,
    sample_out_dir: Path,
    platform: str,
    threads: int,
    config_path: Path | None = None,
    src_dir: Path | None = None,
    report: bool = False,
    skip_existing: bool = False,
) -> int:
    """Run muconespan on a single sample."""
    if skip_existing and (sample_out_dir / "summary.json").exists():
        return 0
    sample_out_dir.mkdir(parents=True, exist_ok=True)
    platform_key = "hifi" if ("hifi" in platform or "pacbio" in platform) else "ont"
    model = MODELS[platform_key]

    cmd = [
        sys.executable,
        "-c",
        "from muc_one_span.cli import main; main()",
    ]
    if config_path:
        cmd.extend(["--config", str(config_path)])
    cmd.extend(
        [
            "run",
            "--input",
            str(fastq_path),
            "--output-dir",
            str(sample_out_dir),
            "--platform",
            platform_key,
            "--clair3-model",
            model,
            "--threads",
            str(threads),
        ]
    )
    if report:
        cmd.extend(["--report", "--report-igv", "embedded"])

    env = dict(os.environ)
    env["PATH"] = f"/home/bernt-popp/miniforge3/envs/env_clair3/bin:{env.get('PATH', '')}"
    active_src = str(src_dir.resolve()) if src_dir else str(ROOT / "src")
    env["PYTHONPATH"] = f"{active_src}:{env.get('PYTHONPATH', '')}"
    env["CLAIR3_MODEL"] = model
    try:
        res = subprocess.run(cmd, cwd=str(ROOT), env=env, capture_output=True, text=True)
        if res.returncode != 0:
            (sample_out_dir / "stderr.log").write_text(res.stderr)
        return res.returncode
    except Exception:
        # Pipeline records failure in run_status.json
        return 1


def process_entry(
    entry: dict[str, Any],
    args: argparse.Namespace,
    rd: Any,
) -> tuple[dict[str, Any], dict[str, Any]]:
    sample_name = f"{entry['design_name']}_{entry['platform']}"
    sample_out = args.output_dir / sample_name
    reads_rel = entry.get("raw_reads_path") or entry.get("reads_file")
    if not reads_rel:
        raise ValueError(f"Missing read path in ledger entry: {entry}")
    reads_path = Path(reads_rel)
    if reads_path.is_absolute() and reads_path.exists():
        raw_reads = reads_path
    elif (args.data_dir / reads_path).exists():
        raw_reads = args.data_dir / reads_path
    else:
        raw_reads = ROOT / reads_path

    truth_rel = entry.get("truth_dir") or (
        Path(entry["truth_fa"]).parent if "truth_fa" in entry else None
    )
    if truth_rel:
        t_path = Path(truth_rel)
        if t_path.is_absolute() and t_path.exists():
            truth_dir = t_path
        elif (args.data_dir / t_path).exists():
            truth_dir = args.data_dir / t_path
        else:
            truth_dir = ROOT / t_path
    else:
        truth_dir = args.data_dir / "truth" / entry["design_name"]

    t0 = time.monotonic()
    exit_code = run_sample(
        sample_name,
        raw_reads,
        sample_out,
        entry["platform"],
        args.threads,
        args.config,
        args.src_dir,
        report=args.report,
        skip_existing=args.skip_existing,
    )
    elapsed = time.monotonic() - t0

    # Offline evaluation against truth
    truth_obj = load_truth(truth_dir, rd)
    obs_obj = load_observation(sample_out, {"exit_code": exit_code, "platform": entry["platform"]})
    eval_row = evaluate_sample(truth_obj, obs_obj)
    eval_row["sample"] = sample_name
    eval_row["wall_seconds"] = round(elapsed, 2)
    eval_row["platform"] = entry["platform"]
    eval_row["category"] = entry["category"]

    # Causal failure diagnosis
    diagnosis = diagnose_first_failure(sample_out, truth_obj, eval_row)
    diagnosis["sample"] = sample_name
    diagnosis["platform"] = entry["platform"]
    diagnosis["category"] = entry["category"]
    diagnosis["lengths"] = entry["lengths"]
    diagnosis["mutation"] = entry["mutation"]
    diagnosis["wall_seconds"] = round(elapsed, 2)

    stage = diagnosis["first_failure_stage"]
    status_sym = "✓" if stage == "none" else f"✗ ({stage})"
    print(f"[{sample_name}] Result: {status_sym} in {elapsed:.1f}s")
    return eval_row, diagnosis


def main() -> int:
    args = parse_args().parse_args()
    if not args.ledger.exists():
        print(f"Error: Ledger file not found: {args.ledger}", file=sys.stderr)
        return 1

    entries = []
    with args.ledger.open() as f:
        for line in f:
            if line.strip():
                entries.append(json.loads(line))

    # Filter by split and platform
    filtered = [e for e in entries if e["split"] == args.split]
    if args.platform != "all":
        filtered = [
            e
            for e in filtered
            if e["platform"] == args.platform
            or (args.platform == "hifi" and ("hifi" in e["platform"] or "pacbio" in e["platform"]))
            or (args.platform == "ont" and "ont" in e["platform"])
        ]
    if args.limit:
        filtered = filtered[: args.limit]

    print(
        f"=== Starting Evaluation for split='{args.split}', platform='{args.platform}' ({len(filtered)} datasets) ==="
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rd = load_repeat_dictionary(_bundled_repeats_path())
    eval_results = []
    failure_atlas = []

    if args.parallel_samples > 1:
        print(f"Running with {args.parallel_samples} parallel workers...")
        with ThreadPoolExecutor(max_workers=args.parallel_samples) as executor:
            futures = [executor.submit(process_entry, entry, args, rd) for entry in filtered]
            for idx, fut in enumerate(futures, 1):
                eval_row, diagnosis = fut.result()
                eval_results.append(eval_row)
                failure_atlas.append(diagnosis)
                print(f"Completed {idx}/{len(filtered)}")
    else:
        for idx, entry in enumerate(filtered, 1):
            print(f"[{idx}/{len(filtered)}] Running {entry['design_name']}_{entry['platform']}...")
            eval_row, diagnosis = process_entry(entry, args, rd)
            eval_results.append(eval_row)
            failure_atlas.append(diagnosis)

    # Sort deterministically
    eval_results.sort(key=lambda r: str(r.get("sample", "")))
    failure_atlas.sort(key=lambda d: str(d.get("sample", "")))

    # Aggregate evaluation
    agg_report = aggregate(eval_results)
    report_file = args.output_dir / "evaluation_report.json"
    report_file.write_text(json.dumps(agg_report, indent=2) + "\n")

    atlas_file = args.output_dir / "failure_atlas.json"
    atlas_file.write_text(json.dumps(failure_atlas, indent=2) + "\n")

    # Build Markdown failure summary
    md_lines = [
        f"# Failure Atlas & Stage Analysis ({args.split.upper()} split)",
        "",
        f"**Total evaluated:** {len(eval_results)}",
        f"**Correct reconstructions:** {sum(1 for d in failure_atlas if d['first_failure_stage'] == 'none')}",
        "",
        "## First Failure Stage Breakdown",
        "",
        "| First Failure Stage | HiFi Count | ONT Count | Total |",
        "|---|---:|---:|---:|",
    ]
    stages = sorted({d["first_failure_stage"] for d in failure_atlas})
    for st in stages:
        h_cnt = sum(
            1
            for d in failure_atlas
            if d["first_failure_stage"] == st
            and ("hifi" in d["platform"] or "pacbio" in d["platform"])
        )
        o_cnt = sum(
            1 for d in failure_atlas if d["first_failure_stage"] == st and "ont" in d["platform"]
        )
        md_lines.append(f"| `{st}` | {h_cnt} | {o_cnt} | {h_cnt + o_cnt} |")

    md_lines.extend(
        [
            "",
            "## Per-Sample Stage Trace",
            "",
            "| Sample | Category | Lengths | Mutation | Stage | Reason | Time (s) |",
            "|---|---|---|---|---|---|---:|",
        ]
    )
    for d in failure_atlas:
        md_lines.append(
            f"| `{d['sample']}` | {d['category']} | {d['lengths']} | {d['mutation']} | `{d['first_failure_stage']}` | {d['reason']} | {d['wall_seconds']} |"
        )

    atlas_md = args.output_dir / "failure_atlas.md"
    atlas_md.write_text("\n".join(md_lines) + "\n")

    print(f"\nCompleted! Evaluation Report: {report_file}")
    print(f"Failure Atlas: {atlas_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
