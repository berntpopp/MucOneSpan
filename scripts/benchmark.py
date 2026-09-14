#!/usr/bin/env python3
"""Benchmark pipeline stage timings on MucOneUp test samples.

Usage:
    export PATH="/path/to/conda/env/bin:$PATH"
    python scripts/benchmark.py [--data-dir tests/data/generated]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark pipeline stages")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("tests/data/generated"),
        help="Directory containing generated test samples",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tests/results/benchmarks"),
        help="Output directory for benchmark results",
    )
    args = parser.parse_args()

    if not args.data_dir.exists():
        print(f"Error: {args.data_dir} not found. Generate test data first.")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)

    from open_pacmuci.config import load_repeat_dictionary
    from open_pacmuci.tools import check_tools

    check_tools(["minimap2", "samtools", "bcftools", "run_clair3.sh"])
    rd = load_repeat_dictionary()

    samples = sorted(args.data_dir.glob("sample_*"))
    results = []

    for sample_dir in samples:
        bam_files = list(sample_dir.glob("*.bam"))
        if not bam_files:
            continue

        bam = bam_files[0]
        sample_name = sample_dir.name
        out = args.output_dir / sample_name
        out.mkdir(parents=True, exist_ok=True)

        timings: dict[str, float] = {}
        print(f"\n{sample_name}:")

        from open_pacmuci.ladder import generate_ladder_fasta
        from open_pacmuci.mapping import get_idxstats, map_reads

        ref = out / "ladder.fa"
        t0 = time.perf_counter()
        generate_ladder_fasta(rd, ref)
        mapped = map_reads(bam, ref, out, threads=4)
        timings["mapping"] = time.perf_counter() - t0
        print(f"  mapping: {timings['mapping']:.2f}s")

        from open_pacmuci.alleles import detect_alleles, parse_idxstats

        t0 = time.perf_counter()
        idxstats = get_idxstats(mapped)
        counts = parse_idxstats(idxstats)
        alleles_result = detect_alleles(counts, min_coverage=10, bam_path=mapped)
        timings["alleles"] = time.perf_counter() - t0
        print(f"  alleles: {timings['alleles']:.2f}s")

        from open_pacmuci.calling import call_variants_per_allele

        t0 = time.perf_counter()
        vcf_paths = call_variants_per_allele(mapped, ref, alleles_result, out, threads=4)
        timings["calling"] = time.perf_counter() - t0
        print(f"  calling: {timings['calling']:.2f}s")

        from open_pacmuci.consensus import build_consensus_per_allele

        t0 = time.perf_counter()
        consensus = build_consensus_per_allele(ref, vcf_paths, alleles_result, out, repeat_dict=rd)
        timings["consensus"] = time.perf_counter() - t0
        print(f"  consensus: {timings['consensus']:.2f}s")

        from open_pacmuci.classify import classify_sequence

        t0 = time.perf_counter()
        for fa_path in consensus.values():
            fa_lines = fa_path.read_text().strip().splitlines()
            sequence = "".join(line for line in fa_lines if not line.startswith(">"))
            classify_sequence(sequence, rd)
        timings["classify"] = time.perf_counter() - t0
        print(f"  classify: {timings['classify']:.2f}s")

        timings["total"] = sum(timings.values())
        print(f"  TOTAL: {timings['total']:.2f}s")

        results.append({"sample": sample_name, "timings": timings})

    summary_path = args.output_dir / "benchmark_results.json"
    summary_path.write_text(json.dumps(results, indent=2) + "\n")
    print(f"\nResults written to {summary_path}")


if __name__ == "__main__":
    main()
