"""Full pipeline execution, separated from Click command declarations."""

from __future__ import annotations

import json
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

from muc_one_span.cli_settings import effective_run_settings, write_run_configuration
from muc_one_span.deprecations import warn_deprecations
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings

if TYPE_CHECKING:
    from muc_one_span.config import RepeatDictionary


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
    engine: str | None = None,
    assay: str | None = None,
) -> None:
    """Run the full MucOneSpan pipeline."""
    from muc_one_span.alleles import detect_alleles, parse_idxstats
    from muc_one_span.calling import call_variants_per_allele
    from muc_one_span.cli import PLATFORM_PRESETS, _bundled_reference
    from muc_one_span.config import load_repeat_dictionary
    from muc_one_span.consensus import build_consensus_per_allele
    from muc_one_span.mapping import get_idxstats, map_reads
    from muc_one_span.pipeline_tail import finish_run
    from muc_one_span.selection_qc import annotate_selection_qc
    from muc_one_span.tools import check_tools, get_tool_versions

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
        **{k: v for k, v in (("engine", engine), ("assay", assay)) if v is not None},
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
    warn_deprecations(settings)
    if settings.run.engine == "hybrid":
        _run_hybrid(out, input_path, rd, settings, configuration_record, report)
        return
    igv_requested = settings.run.report_igv != "off"
    check_tools(
        ["minimap2", "samtools", "bcftools", "run_clair3.sh"]
        + (["create_report"] if igv_requested else [])
    )
    if igv_requested:
        from muc_one_span.report_igv import preflight_igv_report

        preflight_igv_report(settings.run.report_igv, out)
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
        repeat_length_bp=rd.repeat_length_bp,
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

    finish_run(
        out=out,
        input_path=input_path,
        rd=rd,
        settings=settings,
        alleles_result=alleles_result,
        consensus_paths=consensus_paths,
        vcf_paths=vcf_paths,
        tool_versions=tool_versions,
        configuration_record=configuration_record,
        report=report,
        bam_path=bam,
        fasta_path=ref,
    )


def _run_hybrid(
    out: Path,
    input_path: str,
    rd: RepeatDictionary,
    settings: RuntimeSettings,
    configuration_record: dict[str, Any],
    report: bool,
) -> None:
    """Hybrid engine: no mapping, Clair3 or VCF; evidence comes from the reads."""
    from muc_one_span.hybrid.engine import (
        annotate_read_support,
        extra_versions,
        reconstruct_alleles,
    )
    from muc_one_span.pipeline_tail import finish_run
    from muc_one_span.tools import check_tools

    if settings.run.report_igv != "off":
        raise click.BadParameter(
            "IGV tracks are not available for the hybrid engine", param_hint="--report-igv"
        )
    check_tools(["samtools"] if Path(input_path).suffix == ".bam" else [])
    click.echo("Hybrid engine: reconstructing alleles from reads...")
    hybrid = reconstruct_alleles(Path(input_path), out, rd, settings)
    (out / "alleles.json").write_text(json.dumps(hybrid.alleles, indent=2) + "\n")
    finish_run(
        out=out,
        input_path=input_path,
        rd=rd,
        settings=settings,
        alleles_result=hybrid.alleles,
        consensus_paths=hybrid.consensus_paths,
        vcf_paths={},
        tool_versions=extra_versions(hybrid.block["poa_backend"]),
        configuration_record=configuration_record,
        report=report,
        bam_path=None,
        fasta_path=out / "hybrid_references.fa",
        annotate=partial(
            annotate_read_support, rd=rd, members=hybrid.members, settings=settings.hybrid
        ),
        extra_summary={"hybrid": hybrid.block},
    )
