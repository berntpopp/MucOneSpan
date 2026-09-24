"""Full pipeline execution, separated from Click command declarations."""

from __future__ import annotations

import json
from pathlib import Path

import click

from muc_one_span.cli_settings import effective_run_settings, write_run_configuration
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings
from muc_one_span.version import __version__


def execute_pipeline(
    input_path: str,
    output_dir: str,
    reference: str | None,
    clair3_model: str,
    threads: int,
    min_coverage: int,
    min_qual: float,
    report: bool,
    platform: str,
    minimap2_preset: str | None,
    *,
    report_igv: str = "off",
    mapping_timeout: float | None = None,
    settings: RuntimeSettings | None = None,
    configuration: Path | None = None,
) -> None:
    """Run the full MucOneSpan pipeline."""
    from muc_one_span.alleles import detect_alleles, parse_idxstats
    from muc_one_span.calling import call_variants_per_allele
    from muc_one_span.classify import classify_sequence
    from muc_one_span.cli import PLATFORM_PRESETS, _bundled_reference
    from muc_one_span.config import load_repeat_dictionary
    from muc_one_span.consensus import build_consensus_per_allele
    from muc_one_span.mapping import get_idxstats, map_reads
    from muc_one_span.selection_qc import annotate_selection_qc
    from muc_one_span.tools import check_tools, get_tool_versions
    from muc_one_span.vcf import parse_vcf_variants

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "run_configuration.json").unlink(missing_ok=True)
    # Validate files within the recorded invocation so a failed rerun cannot
    # leave a stale completed sidecar. Preserve Click's usage-error exit code.
    for option, value in (("--input", input_path), ("--reference", reference)):
        if value is not None and not Path(value).is_file():
            raise click.BadParameter(f"File does not exist: {value}", param_hint=option)
    preset = minimap2_preset or PLATFORM_PRESETS[platform]

    ref = Path(reference) if reference else _bundled_reference()
    settings = effective_run_settings(
        settings or DEFAULT_SETTINGS,
        reference=reference,
        clair3_model=clair3_model,
        threads=threads,
        min_coverage=min_coverage,
        min_qual=min_qual,
        report=report,
        report_igv=report_igv,
        platform=platform,
        minimap2_preset=minimap2_preset,
        mapping_timeout=mapping_timeout
        if mapping_timeout is not None
        else (settings or DEFAULT_SETTINGS).run.mapping_timeout,
    )
    if reference is None and (
        settings.repeat_dictionary is not None
        or settings.reference_layout.pre != DEFAULT_SETTINGS.reference_layout.pre
        or settings.reference_layout.after != DEFAULT_SETTINGS.reference_layout.after
        or settings.consensus.flank_length != DEFAULT_SETTINGS.consensus.flank_length
    ):
        raise click.BadParameter(
            "A configured dictionary, reference layout or flank length requires an explicit matching reference.",
            param_hint="--reference",
        )
    rd = load_repeat_dictionary(
        Path(settings.repeat_dictionary) if settings.repeat_dictionary else None
    )
    try:
        settings.reference_layout.validate_repeats(rd.repeats)
        settings.consensus.validate_flanks(rd.flanking_left, rd.flanking_right)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--config") from exc
    configuration_record = write_run_configuration(
        settings, configuration, Path(input_path), ref, out
    )
    check_tools(["minimap2", "samtools", "bcftools", "run_clair3.sh"])
    tool_versions = get_tool_versions(["minimap2", "samtools", "bcftools", "run_clair3.sh"])

    # Step 1: Map reads
    click.echo("Step 1/5: Mapping reads...")
    bam = map_reads(
        Path(input_path), ref, out, threads, preset=preset, timeout=settings.run.mapping_timeout
    )

    # Step 2: Detect alleles
    click.echo("Step 2/5: Detecting alleles...")
    idxstats = get_idxstats(bam)
    counts = parse_idxstats(idxstats)
    alleles_result = detect_alleles(
        counts,
        min_coverage,
        bam_path=bam,
        settings=settings.allele_selection,
        reference_layout=settings.reference_layout,
        platform=settings.run.platform,
    )
    annotate_selection_qc(alleles_result, settings.allele_selection)
    (out / "alleles.json").write_text(json.dumps(alleles_result, indent=2) + "\n")
    click.echo(f"  Alleles: {alleles_result}")

    # Step 3: Call variants
    click.echo("Step 3/5: Calling variants...")
    vcf_paths = call_variants_per_allele(
        bam,
        ref,
        alleles_result,
        out,
        clair3_model,
        threads,
        min_qual=min_qual,
        platform=platform,
        preset=preset,
        settings=settings.calling,
        read_phasing_settings=settings.read_phasing,
    )

    # Step 4: Build consensus
    click.echo("Step 4/5: Building consensus...")
    consensus_paths = build_consensus_per_allele(
        ref,
        vcf_paths,
        alleles_result,
        out,
        repeat_dict=rd,
        settings=settings.consensus,
        reference_layout=settings.reference_layout,
    )

    # Step 5: Classify repeats
    click.echo("Step 5/5: Classifying repeats...")
    from muc_one_span.classify import validate_mutations_against_vcf

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
            bam_path=bam,
            vcf_paths=vcf_paths,
            fasta_path=ref,
            execution_status={"status": "analysis_completed"},
        )
        click.echo(f"Report: {report_path}")

    summary["run_status"] = {"status": "completed"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    click.echo("Pipeline complete.")
