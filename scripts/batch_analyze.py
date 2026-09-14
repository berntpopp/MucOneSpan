#!/usr/bin/env python3
"""Batch-run open-pacmuci pipeline on all generated test samples and report results.

Produces a TP/FP/FN analysis comparing detected mutations against ground truth
from MucOneUp's vntr_structure.txt files.

Usage:
    export PATH="/home/bernt-popp/miniforge3/envs/env_pacbio/bin:$PATH"
    export CLAIR3_MODEL=/path/to/clair3/models/hifi
    python scripts/batch_analyze.py [samples_dir] [output_dir]
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

SAMPLES_DIR = Path("tests/data/generated")
OUTPUT_DIR = Path("tests/results")
CLAIR3_MODEL = os.environ.get("CLAIR3_MODEL", "")


def parse_ground_truth(sample_dir: Path) -> dict:
    """Parse MucOneUp vntr_structure.txt to extract expected mutations."""
    structure_files = sorted(sample_dir.glob("*.vntr_structure.txt"))
    if not structure_files:
        return {"mutation": "unknown", "haplotypes": {}}

    text = structure_files[0].read_text()
    # Parse mutation name from header
    mut_match = re.search(r"Mutation Applied:\s+(\S+)", text)
    mutation = mut_match.group(1) if mut_match else "normal"

    haplotypes = {}
    for line in text.strip().splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            name = parts[0].strip()
            structure = parts[1].strip()
            units = [u for u in structure.split("-") if u]
            has_mutation = any(u.endswith("m") for u in units)  # "Xm" = mutated repeat
            n_repeats = len(units)
            haplotypes[name] = {
                "structure": structure,
                "has_mutation": has_mutation,
                "n_repeats": n_repeats,
            }

    return {"mutation": mutation, "haplotypes": haplotypes}


def run_pipeline(sample_dir: Path, output_dir: Path, platform: str = "hifi") -> dict | None:
    """Run open-pacmuci pipeline on a sample, return parsed results."""
    # Find reads BAM
    bams = sorted(sample_dir.glob("*_reads_amplicon_aligned.bam"))
    if not bams:
        # Try FASTQ
        fastqs = sorted(sample_dir.glob("*_reads*.fastq")) + sorted(sample_dir.glob("*_reads*.fq"))
        if not fastqs:
            return None
        input_file = fastqs[0]
    else:
        input_file = bams[0]

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "uv",
        "run",
        "open-pacmuci",
        "run",
        "--input",
        str(input_file),
        "--output-dir",
        str(output_dir),
        "--threads",
        "4",
        "--min-qual",
        "5.0",
        "--clair3-model",
        CLAIR3_MODEL,
        "--platform",
        platform,
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )
        if result.returncode != 0:
            return {"error": result.stderr[-500:] if result.stderr else "unknown error"}
    except subprocess.TimeoutExpired:
        return {"error": "timeout"}

    # Parse results
    summary_file = output_dir / "summary.json"
    if not summary_file.exists():
        return {"error": "no summary.json produced"}

    summary: dict = json.loads(summary_file.read_text())
    return summary


def analyze_results(
    sample_name: str,
    ground_truth: dict,
    pipeline_result: dict | None,
) -> dict:
    """Compare pipeline results against ground truth."""
    analysis = {
        "sample": sample_name,
        "expected_mutation": ground_truth["mutation"],
        "status": "error",
        "details": {},
    }

    if pipeline_result is None:
        analysis["details"]["error"] = "no reads found"
        return analysis

    if "error" in pipeline_result:
        analysis["details"]["error"] = pipeline_result["error"]
        return analysis

    # Check allele detection
    alleles = pipeline_result.get("alleles", {})
    analysis["details"]["alleles_detected"] = bool(alleles.get("allele_1"))

    # Check classifications
    classifications = pipeline_result.get("classifications", {})
    detected_mutations = []

    for allele_key, classif in classifications.items():
        mutations = classif.get("mutations", [])
        for mut in mutations:
            mut_info = {
                "allele": allele_key,
                "repeat_index": mut.get("repeat_index"),
                "template_match": mut.get("template_match", False),
                "mutation_name": mut.get("mutation_name", "unknown"),
                "vcf_support": mut.get("vcf_support"),
                "vcf_qual": mut.get("vcf_qual"),
                "boundary": mut.get("boundary", False),
            }
            detected_mutations.append(mut_info)

    analysis["details"]["detected_mutations"] = detected_mutations

    # Classify as TP/FP/FN
    expected_mut = ground_truth["mutation"]
    has_expected = expected_mut != "normal"

    # Check if any detected mutation matches expected
    template_matches = [
        m
        for m in detected_mutations
        if m.get("template_match") and expected_mut in str(m.get("mutation_name", ""))
    ]
    any_mutation_detected = len(detected_mutations) > 0

    if has_expected:
        if template_matches:
            analysis["status"] = "TP"
        elif any_mutation_detected:
            # Detected something but not the right template
            analysis["status"] = "TP_partial"
        else:
            analysis["status"] = "FN"
    else:
        if any_mutation_detected:
            analysis["status"] = "FP"
        else:
            analysis["status"] = "TN"

    # Add confidence info
    for allele_key, classif in classifications.items():
        structure = classif.get("structure", "")
        analysis["details"][f"{allele_key}_structure_len"] = len(structure.split())

    return analysis


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Batch analyze test samples")
    parser.add_argument("samples_dir", nargs="?", default=str(SAMPLES_DIR))
    parser.add_argument("output_dir", nargs="?", default=str(OUTPUT_DIR))
    parser.add_argument(
        "--platform",
        default="hifi",
        choices=["hifi", "ont"],
        help="Sequencing platform (default: hifi)",
    )
    args = parser.parse_args()

    samples_dir = Path(args.samples_dir)
    output_base = Path(args.output_dir)

    # Get all sample directories that have reads
    sample_dirs = sorted(
        d
        for d in samples_dir.iterdir()
        if d.is_dir()
        and (list(d.glob("*_reads_amplicon_aligned.bam")) or list(d.glob("*_reads*.fastq")))
    )

    print(f"Found {len(sample_dirs)} samples to analyze")
    print(f"Output: {output_base}")
    print()

    results = []
    for i, sample_dir in enumerate(sample_dirs, 1):
        name = sample_dir.name
        print(f"[{i}/{len(sample_dirs)}] {name}...", end=" ", flush=True)

        ground_truth = parse_ground_truth(sample_dir)
        out_dir = output_base / name

        pipeline_result = run_pipeline(sample_dir, out_dir, platform=args.platform)
        analysis = analyze_results(name, ground_truth, pipeline_result)
        results.append(analysis)

        print(f"{analysis['status']} (expected: {ground_truth['mutation']})")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    tp = sum(1 for r in results if r["status"] == "TP")
    tp_partial = sum(1 for r in results if r["status"] == "TP_partial")
    fn = sum(1 for r in results if r["status"] == "FN")
    fp = sum(1 for r in results if r["status"] == "FP")
    tn = sum(1 for r in results if r["status"] == "TN")
    errors = sum(1 for r in results if r["status"] == "error")

    print(f"  TP (template match):  {tp}")
    print(f"  TP (partial):         {tp_partial}")
    print(f"  FN (missed):          {fn}")
    print(f"  FP (false alarm):     {fp}")
    print(f"  TN (correct normal):  {tn}")
    print(f"  Errors:               {errors}")
    print()

    if fn > 0:
        print("FALSE NEGATIVES:")
        for r in results:
            if r["status"] == "FN":
                print(f"  {r['sample']}: expected {r['expected_mutation']}")
        print()

    if fp > 0:
        print("FALSE POSITIVES:")
        for r in results:
            if r["status"] == "FP":
                det = r["details"].get("detected_mutations", [])
                print(f"  {r['sample']}: detected {len(det)} mutations")
                for m in det:
                    print(
                        f"    repeat {m['repeat_index']}: "
                        f"vcf_qual={m.get('vcf_qual')}, "
                        f"boundary={m.get('boundary')}"
                    )
        print()

    if errors > 0:
        print("ERRORS:")
        for r in results:
            if r["status"] == "error":
                err = r["details"].get("error", "unknown")
                print(f"  {r['sample']}: {err[:200]}")
        print()

    # Write full results JSON
    output_base.mkdir(parents=True, exist_ok=True)
    results_file = output_base / "batch_results.json"
    results_file.write_text(json.dumps(results, indent=2, default=str) + "\n")
    print(f"Full results: {results_file}")


if __name__ == "__main__":
    main()
