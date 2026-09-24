"""Standalone report input loading and explicit execution provenance."""

from __future__ import annotations

import json
from pathlib import Path

import click


def execute_report(
    input_path: str,
    output: str,
    sample_name: str | None,
    repeats: str | None,
    report_igv: str = "off",
    fasta: str | None = None,
    bam: str | None = None,
    vcf: str | None = None,
    allele_vcf: tuple[tuple[str, str], ...] = (),
    run_status: str | None = None,
) -> None:
    """Render a summary using its own sidecar, never the report destination's."""
    try:
        from muc_one_span.report import generate_report
    except ImportError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e

    try:
        summary = json.loads(Path(input_path).read_text())
    except (json.JSONDecodeError, ValueError) as e:
        click.echo(f"Error: Failed to parse JSON from {input_path}: {e}", err=True)
        raise SystemExit(1) from e

    name = sample_name or Path(input_path).parent.name

    detailed = None
    if repeats:
        try:
            detailed = json.loads(Path(repeats).read_text())
        except (json.JSONDecodeError, ValueError) as e:
            click.echo(f"Error: Failed to parse JSON from {repeats}: {e}", err=True)
            raise SystemExit(1) from e

    status_path = Path(run_status) if run_status else Path(input_path).parent / "run_status.json"
    execution_status = None
    if run_status or status_path.exists():
        try:
            execution_status = json.loads(status_path.read_text())
            if not isinstance(execution_status, dict):
                raise ValueError("status must be a JSON object")
        except (OSError, ValueError) as error:
            raise click.ClickException(
                f"Failed to read execution status {status_path}: {error}"
            ) from error

    vcf_paths: dict[str, Path] | None = None
    if allele_vcf:
        vcf_paths = {}
        for label, supplied in allele_vcf:
            if label in vcf_paths:
                raise click.BadParameter(f"Duplicate VCF label: {label}", param_hint="--allele-vcf")
            vcf_paths[label] = Path(supplied)
    elif not vcf:
        # Consensus records an absolute VCF path even when calling received a
        # relative output directory. Prefer this reproducible provenance. Legacy
        # relative paths without context are interpreted beside the summary.
        for key, info in summary.get("alleles", {}).items():
            if not isinstance(info, dict):
                continue
            context = info.get("consensus_context") or {}
            recorded = context.get("vcf_path") or info.get("vcf_path")
            if recorded:
                path = Path(recorded)
                if not path.is_absolute():
                    path = Path(input_path).parent / path
                if vcf_paths is None:
                    vcf_paths = {}
                vcf_paths[key] = path

    try:
        out_path = generate_report(
            summary,
            Path(output),
            sample_name=name,
            detailed_repeats=detailed,
            report_igv=report_igv,
            fasta_path=Path(fasta) if fasta else None,
            bam_path=Path(bam) if bam else None,
            vcf_path=Path(vcf) if vcf else None,
            vcf_paths=vcf_paths,
            execution_status=execution_status,
        )
    except ValueError as error:
        # A bad or unknown recorded clinical_decision section
        # (decision_settings.resolve_decision_settings) must fail cleanly.
        raise click.ClickException(str(error)) from error
    click.echo(f"Report written to {out_path}")
