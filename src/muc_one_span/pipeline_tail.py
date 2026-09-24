"""Shared classification, summary and report tail for both reconstruction engines."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import click

from muc_one_span.config import RepeatDictionary
from muc_one_span.settings import RuntimeSettings
from muc_one_span.version import __version__

Annotate = Callable[[str, str, dict[str, Any]], dict[str, Any]]


def finish_run(
    *,
    out: Path,
    input_path: str,
    rd: RepeatDictionary,
    settings: RuntimeSettings,
    alleles_result: dict[str, Any],
    consensus_paths: dict[str, Path],
    vcf_paths: dict[str, Path],
    tool_versions: dict[str, str],
    configuration_record: dict[str, Any],
    report: bool,
    bam_path: Path | None,
    fasta_path: Path | None,
    annotate: Annotate | None = None,
    extra_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Classify each consensus, write repeats/summary files and the optional report."""
    # Imported here so that tests patching these module attributes keep working.
    from muc_one_span.classify import classify_sequence, validate_mutations_against_vcf
    from muc_one_span.vcf import parse_vcf_variants

    # Step 5: Classify repeats
    click.echo("Step 5/5: Classifying repeats...")
    all_results: dict[str, dict] = {}
    for allele_key, fa_path in consensus_paths.items():
        fa_lines = fa_path.read_text().strip().splitlines()
        sequence = "".join(line for line in fa_lines if not line.startswith(">"))
        result = classify_sequence(sequence, rd, settings=settings.classification)

        # VCF-backed validation if VCF available
        if allele_key in vcf_paths:
            vcf_variants = parse_vcf_variants(vcf_paths[allele_key])
            result = validate_mutations_against_vcf(
                result,
                vcf_variants=vcf_variants,
                sequence=sequence,
                repeat_dict=rd,
                consensus_context=alleles_result[allele_key].get("consensus_context"),
                settings=settings.confidence,
            )
        if annotate is not None:
            result = annotate(allele_key, sequence, result)

        all_results[allele_key] = result
        click.echo(f"  {allele_key}: {result['structure']}")
        if result.get("allele_confidence") is not None:
            click.echo(f"    confidence: {result['allele_confidence']:.2f}")

    (out / "alleles.json").write_text(json.dumps(alleles_result, indent=2) + "\n")

    # Write combined outputs
    (out / "repeats.json").write_text(json.dumps(all_results, indent=2) + "\n")
    structures = {k: v["structure"] for k, v in all_results.items()}
    (out / "repeats.txt").write_text("\n".join(f"{k}: {v}" for k, v in structures.items()) + "\n")

    # Summary
    summary = {
        "run_status": {"status": "analysis_completed"},
        "alleles": alleles_result,
        "classifications": {
            k: {
                "structure": v["structure"],
                "mutations": v["mutations_detected"],
                "reconstruction_status": v.get("reconstruction_status", "unverified"),
                "ambiguous_bases": v.get("ambiguous_bases", 0),
                "classification_coverage": v.get("classification_coverage"),
                "vcf_projection": v.get("vcf_projection"),
                "confidence_semantics": "heuristic_dictionary_fit_not_probability",
            }
            for k, v in all_results.items()
        },
        "tool_versions": tool_versions,
        "pipeline_version": __version__,
        "configuration": configuration_record,
        **(extra_summary or {}),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    effective_igv = settings.run.report_igv
    if report or effective_igv != "off":
        from muc_one_span.report import generate_report

        report_path = out / "report.html"
        generate_report(
            summary,
            report_path,
            sample_name=Path(input_path).stem,
            detailed_repeats=all_results,
            report_igv=effective_igv,
            bam_path=bam_path,
            vcf_paths=vcf_paths,
            fasta_path=fasta_path,
            execution_status={"status": "analysis_completed"},
        )
        click.echo(f"Report: {report_path}")

    summary["run_status"] = {"status": "completed"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    click.echo("Pipeline complete.")
    return summary
