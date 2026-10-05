"""Full pipeline execution, separated from Click command declarations."""

from __future__ import annotations

import json
import shutil
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

from muc_one_span.cli_settings import effective_run_settings, write_run_configuration
from muc_one_span.deprecations import (
    LADDER_ENGINE,
    explicit_command_line_options,
    ignored_options,
    warn_deprecations,
    warn_ignored,
)
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings

if TYPE_CHECKING:
    from muc_one_span.config import RepeatDictionary
    from muc_one_span.pipeline_tail import IgvInputs

# Hybrid --report-igv: display-only alignment of the reads to the polished consensus.
HYBRID_IGV_TOOLS = ["minimap2", "samtools"]
HYBRID_IGV_DIR = "igv"


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
    igv_session: bool | None = None,
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
        **{
            k: v
            for k, v in (("engine", engine), ("assay", assay), ("igv_session", igv_session))
            if v is not None
        },
    )
    # Only the ladder aligns to a reference FASTA; the hybrid engine builds its own.
    if (
        reference is None
        and settings.run.engine == LADDER_ENGINE
        and (
            settings.repeat_dictionary is not None
            or settings.reference_layout.pre != DEFAULT_SETTINGS.reference_layout.pre
            or settings.reference_layout.after != DEFAULT_SETTINGS.reference_layout.after
            or settings.consensus.flank_length != DEFAULT_SETTINGS.consensus.flank_length
        )
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
    if settings.run.engine == LADDER_ENGINE and settings.run.igv_session:
        raise click.BadParameter(
            "IGV session output is available for the hybrid engine only; the ladder "
            "engine already writes mapping.bam",
            param_hint="--igv-session",
        )
    ignored = ignored_options(settings, explicit_command_line_options())
    warn_deprecations(settings)
    warn_ignored(ignored)
    configuration_record = write_run_configuration(
        settings, configuration, Path(input_path), ref, out, ignored_options=ignored
    )
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
    """Hybrid engine: no Clair3 or VCF; evidence comes from the reads.

    With ``--report-igv`` the reads are additionally aligned (display only) to each
    allele with its mutated units restored to their canonical parents, so reads
    carrying a detected mutation show it, and the mutations become an IGV track.
    """
    from muc_one_span.hybrid.engine import (
        annotate_read_support,
        extra_versions,
        reconstruct_alleles,
    )
    from muc_one_span.pipeline_tail import finish_run
    from muc_one_span.tools import check_tools, get_tool_versions

    report_igv = settings.run.report_igv != "off"
    igv = report_igv or settings.run.igv_session
    tools = ["samtools"] if Path(input_path).suffix == ".bam" else []
    if igv:
        tools = list(dict.fromkeys([*tools, *HYBRID_IGV_TOOLS]))
    if report_igv:
        tools.append("create_report")
    check_tools(tools)
    if report_igv:
        from muc_one_span.report_igv import preflight_igv_report

        preflight_igv_report(settings.run.report_igv, out)
    click.echo("Hybrid engine: reconstructing alleles from reads...")
    hybrid = reconstruct_alleles(Path(input_path), out, rd, settings)
    versions = extra_versions(hybrid.block["poa_backend"])
    igv_tracks = None
    if igv:
        consensus = {
            allele: "".join(
                line for line in path.read_text().splitlines() if not line.startswith(">")
            )
            for allele, path in hybrid.consensus_paths.items()
        }
        igv_tracks = partial(
            _hybrid_igv_inputs, out, input_path, rd, settings, consensus, hybrid.read_names
        )
        versions.update(get_tool_versions(HYBRID_IGV_TOOLS))
    # alleles.json is written once, by finish_run, after read-support annotation.
    finish_run(
        out=out,
        input_path=input_path,
        rd=rd,
        settings=settings,
        alleles_result=hybrid.alleles,
        consensus_paths=hybrid.consensus_paths,
        vcf_paths={},
        tool_versions=versions,
        configuration_record=configuration_record,
        report=report,
        bam_path=None,
        fasta_path=out / "hybrid_references.fa",
        annotate=partial(
            annotate_read_support, rd=rd, members=hybrid.members, settings=settings.hybrid
        ),
        extra_summary={"hybrid": hybrid.block},
        igv_tracks=igv_tracks,
    )
    if igv and not settings.run.igv_session:
        # The HTML report embeds what it needs; keep igv/ only when it was requested.
        shutil.rmtree(out / HYBRID_IGV_DIR, ignore_errors=True)


def _hybrid_igv_inputs(
    out: Path,
    input_path: str,
    rd: RepeatDictionary,
    settings: RuntimeSettings,
    consensus: dict[str, str],
    read_names: dict[str, list[str]],
    classifications: dict[str, dict[str, Any]],
) -> IgvInputs:
    """Write the IGV display reference and tracks, then align each allele's reads.

    Called after classification: the display reference restores each mutated unit
    to its canonical parent so that reads carrying the mutation show it, and each
    allele shows only the reads the engine assigned to it (at most
    ``hybrid.igv_max_reads_per_allele``, seeded subsample). Files go to ``igv/``
    (``igv_reference.fa`` + ``.fai``, ``loci.bed``, track BEDs, ``mapping.bam`` +
    ``.bai`` and ``igv_session.xml`` for IGV Desktop).
    """
    import random

    from muc_one_span.hybrid.igv_alignment import align_assigned_reads
    from muc_one_span.hybrid.igv_gene import display_flanks, load_gene_model, parse_region
    from muc_one_span.hybrid.igv_session import locus_string, write_igv_session
    from muc_one_span.hybrid.igv_tracks import write_igv_inputs
    from muc_one_span.pipeline_tail import IgvInputs as Inputs

    igv_dir = out / HYBRID_IGV_DIR
    h = settings.hybrid
    gene = load_gene_model(Path(h.igv_gene_annotation) if h.igv_gene_annotation else None)
    available = (len(rd.flanking_left), len(rd.flanking_right))
    flanks = display_flanks(gene, parse_region(rd.vntr_region), h.igv_gene_margin_bp, available)
    files = write_igv_inputs(igv_dir, consensus, classifications, rd, flanks, gene)
    rng = random.Random(h.seed)
    shown = {
        allele: names
        if len(names) <= h.igv_max_reads_per_allele
        else sorted(rng.sample(sorted(names), h.igv_max_reads_per_allele))
        for allele, names in sorted(read_names.items())
    }
    click.echo("Aligning each allele's assigned reads for the IGV report...")
    bam = align_assigned_reads(
        Path(input_path),
        files.fasta,
        shown,
        igv_dir,
        threads=settings.run.threads,
        preset=h.igv_minimap2_preset,
        timeout=settings.run.mapping_timeout,
    )
    window_bp = h.igv_context_units * rd.repeat_length_bp
    rows = [line.split("\t") for line in files.loci.read_text().splitlines()]
    first = next((r for r in rows if r[3] in files.mutation_names), rows[0])
    write_igv_session(
        igv_dir,
        files.fasta,
        files.tracks,
        bam,
        files.loci,
        locus_string(first[0], int(first[1]), int(first[2]), window_bp),
    )
    return Inputs(
        fasta=files.fasta,
        bam=bam,
        loci_bed=files.loci,
        tracks=files.tracks,
        # igv-reports splits --flanking/--window evenly between the two sides: twice
        # the longest contig embeds every contig whole around any row.
        flanking=2 * max(files.contig_lengths.values()),
        window=2 * window_bp,
        preferred_rows=files.mutation_names,
    )
