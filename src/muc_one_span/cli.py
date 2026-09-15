"""Click CLI for MucOneSpan pipeline."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import click

from muc_one_span.cli_settings import (
    configure_context,
    current_configuration_path,
    current_settings,
    validate_stage_options,
)
from muc_one_span.mapping import PLATFORM_PRESETS
from muc_one_span.run_status import record_run_status
from muc_one_span.settings import DEFAULT_SETTINGS
from muc_one_span.version import __version__


@click.group()
@click.version_option(version=__version__, prog_name="muconespan")
@click.option(
    "-v", "--verbose", count=True, help="Increase verbosity (-v for INFO, -vv for DEBUG)."
)
@click.option("-q", "--quiet", is_flag=True, help="Suppress non-error output.")
@click.option(
    "--config",
    "configuration",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Runtime settings JSON; explicit command options override file values.",
)
@click.pass_context
def main(ctx: click.Context, verbose: int, quiet: bool, configuration: Path | None) -> None:
    """MucOneSpan: MUC1 VNTR analysis pipeline for PacBio HiFi and ONT amplicon data."""
    configure_context(ctx, configuration)
    if quiet:
        level = logging.ERROR
    elif verbose >= 2:
        level = logging.DEBUG
    elif verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING

    logging.basicConfig(
        level=level,
        format="%(name)s %(levelname)s: %(message)s",
    )


@main.command()
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    default="reference_ladder.fa",
    help="Output FASTA path.",
)
@click.option(
    "--min-units",
    type=int,
    default=DEFAULT_SETTINGS.reference_layout.min_units,
    help="Minimum repeat units.",
)
@click.option(
    "--max-units",
    type=int,
    default=DEFAULT_SETTINGS.reference_layout.max_units,
    help="Maximum repeat units.",
)
@click.option(
    "--flank-length",
    type=int,
    default=DEFAULT_SETTINGS.consensus.flank_length,
    help="Flanking sequence length (bp).",
)
@click.option(
    "--repeats-db",
    type=click.Path(exists=True),
    default=None,
    help="Custom repeat dictionary JSON.",
)
def ladder(
    output: str,
    min_units: int,
    max_units: int,
    flank_length: int,
    repeats_db: str | None,
) -> None:
    """Generate or regenerate the reference ladder FASTA."""
    from muc_one_span.config import load_repeat_dictionary
    from muc_one_span.ladder import generate_ladder_fasta

    validate_stage_options()
    rd = load_repeat_dictionary(Path(repeats_db) if repeats_db else None)
    settings = current_settings()
    out_path = generate_ladder_fasta(
        rd,
        Path(output),
        min_units,
        max_units,
        flank_length,
        settings=settings.consensus,
        reference_layout=settings.reference_layout,
    )
    click.echo(f"Ladder written to {out_path} ({max_units - min_units + 1} contigs)")


@main.command(name="map")
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ or BAM file.",
)
@click.option(
    "--reference",
    "-r",
    type=click.Path(exists=True),
    default=None,
    help="Reference FASTA (defaults to bundled ladder).",
)
@click.option("--output-dir", "-o", type=click.Path(), default=".", help="Output directory.")
@click.option(
    "--threads", "-t", type=int, default=DEFAULT_SETTINGS.run.threads, help="Number of threads."
)
@click.option(
    "--platform",
    type=click.Choice(["hifi", "ont"], case_sensitive=False),
    default=DEFAULT_SETTINGS.run.platform,
    help="Sequencing platform (default: hifi).",
)
@click.option(
    "--mapping-timeout",
    type=float,
    default=DEFAULT_SETTINGS.run.mapping_timeout,
    help="Total mapping timeout in seconds (finite and positive; default 3600).",
)
@click.option(
    "--minimap2-preset",
    type=str,
    default=None,
    help="minimap2 -x preset (auto-selected from --platform if not set).",
)
def map_cmd(
    input_path: str,
    reference: str | None,
    output_dir: str,
    threads: int,
    platform: str,
    minimap2_preset: str | None,
    mapping_timeout: float = DEFAULT_SETTINGS.run.mapping_timeout,
) -> None:
    """Map reads to the ladder reference with minimap2."""
    from muc_one_span.mapping import map_reads
    from muc_one_span.tools import check_tools

    validate_stage_options()
    check_tools(["minimap2", "samtools"])

    preset = minimap2_preset or PLATFORM_PRESETS[platform]
    ref = Path(reference) if reference else _bundled_reference()
    bam = map_reads(
        Path(input_path), ref, Path(output_dir), threads, preset=preset, timeout=mapping_timeout
    )
    click.echo(f"Mapping written to {bam}")


@main.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Input BAM file (mapped to ladder).",
)
@click.option(
    "--min-coverage",
    type=int,
    default=DEFAULT_SETTINGS.run.min_coverage,
    help="Minimum read coverage.",
)
@click.option("--output-dir", "-o", type=click.Path(), default=".", help="Output directory.")
def alleles(input_path: str, min_coverage: int, output_dir: str) -> None:
    """Determine allele lengths from mapping."""
    from muc_one_span.alleles import detect_alleles, parse_idxstats
    from muc_one_span.mapping import get_idxstats

    validate_stage_options()
    bam = Path(input_path)
    idxstats_output = get_idxstats(bam)
    counts = parse_idxstats(idxstats_output)
    settings = current_settings()
    result = detect_alleles(
        counts,
        min_coverage,
        bam_path=bam,
        settings=settings.allele_selection,
        reference_layout=settings.reference_layout,
        platform=settings.run.platform,
    )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "alleles.json"
    out_file.write_text(json.dumps(result, indent=2) + "\n")
    click.echo(f"Alleles: {result}")


@main.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Input BAM file.",
)
@click.option(
    "--reference",
    "-r",
    required=True,
    type=click.Path(exists=True),
    help="Reference FASTA.",
)
@click.option(
    "--alleles-json",
    "-a",
    required=True,
    type=click.Path(exists=True),
    help="Alleles JSON from 'alleles' command.",
)
@click.option("--output-dir", "-o", type=click.Path(), default=".", help="Output directory.")
@click.option(
    "--clair3-model",
    type=str,
    default=DEFAULT_SETTINGS.run.clair3_model,
    help="Path to Clair3 model.",
)
@click.option(
    "--threads", "-t", type=int, default=DEFAULT_SETTINGS.run.threads, help="Number of threads."
)
@click.option(
    "--min-qual",
    type=float,
    default=DEFAULT_SETTINGS.run.min_qual,
    help="Minimum QUAL score for VCF filtering (default 5.0).",
)
@click.option(
    "--platform",
    type=click.Choice(["hifi", "ont"], case_sensitive=False),
    default=DEFAULT_SETTINGS.run.platform,
    help="Sequencing platform (default: hifi).",
)
@click.option(
    "--minimap2-preset",
    type=str,
    default=None,
    help="minimap2 -x preset (auto-selected from --platform if not set).",
)
def call(
    input_path: str,
    reference: str,
    alleles_json: str,
    output_dir: str,
    clair3_model: str,
    threads: int,
    min_qual: float,
    platform: str,
    minimap2_preset: str | None,
) -> None:
    """Call variants with Clair3."""
    from muc_one_span.calling import call_variants_per_allele
    from muc_one_span.tools import check_tools

    validate_stage_options()
    check_tools(["minimap2", "samtools", "bcftools", "run_clair3.sh"])

    try:
        alleles_data = json.loads(Path(alleles_json).read_text())
    except (json.JSONDecodeError, ValueError) as exc:
        click.echo(f"Error: failed to parse alleles JSON: {exc}", err=True)
        sys.exit(1)
    preset = minimap2_preset or PLATFORM_PRESETS[platform]
    vcfs = call_variants_per_allele(
        Path(input_path),
        Path(reference),
        alleles_data,
        Path(output_dir),
        clair3_model,
        threads,
        min_qual=min_qual,
        platform=platform,
        preset=preset,
        settings=current_settings().calling,
        read_phasing_settings=current_settings().read_phasing,
    )
    (Path(output_dir) / "alleles.json").write_text(json.dumps(alleles_data, indent=2) + "\n")
    for key, vcf in vcfs.items():
        click.echo(f"{key}: {vcf}")


@main.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Input BAM file.",
)
@click.option(
    "--reference",
    "-r",
    required=True,
    type=click.Path(exists=True),
    help="Reference FASTA.",
)
@click.option(
    "--alleles-json",
    "-a",
    required=True,
    type=click.Path(exists=True),
    help="Alleles JSON.",
)
@click.option("--output-dir", "-o", type=click.Path(), default=".", help="Output directory.")
@click.option(
    "--repeats-db",
    type=click.Path(exists=True),
    default=None,
    help="Custom repeat dictionary for anchor-aware trimming.",
)
def consensus(
    input_path: str,
    reference: str,
    alleles_json: str,
    output_dir: str,
    repeats_db: str | None,
) -> None:
    """Build per-allele consensus sequences."""
    from muc_one_span.config import load_repeat_dictionary
    from muc_one_span.consensus import build_consensus_per_allele
    from muc_one_span.tools import check_tools

    rd = load_repeat_dictionary(Path(repeats_db) if repeats_db else None)
    check_tools(["samtools", "bcftools"])

    try:
        alleles_data = json.loads(Path(alleles_json).read_text())
    except (json.JSONDecodeError, ValueError) as exc:
        click.echo(f"Error: failed to parse alleles JSON: {exc}", err=True)
        sys.exit(1)
    out = Path(output_dir)
    vcf_paths: dict[str, Path] = {}
    for key in ["allele_1", "allele_2"]:
        info = alleles_data.get(key, {})
        if info.get("candidate_duplicate_of"):
            continue
        vcf = Path(info["vcf_path"]) if info.get("vcf_path") else out / key / "variants.vcf.gz"
        if vcf.exists():
            vcf_paths[key] = vcf

    if not vcf_paths and (out / "merged" / "variants.vcf.gz").exists():
        vcf_paths["allele_1"] = out / "merged" / "variants.vcf.gz"
    if not vcf_paths:
        click.echo(
            "Warning: no VCF files found. Expected allele_1/variants.vcf.gz "
            "or allele_2/variants.vcf.gz in the output directory.",
            err=True,
        )

    fastas = build_consensus_per_allele(
        Path(reference),
        vcf_paths,
        alleles_data,
        out,
        repeat_dict=rd,
        settings=current_settings().consensus,
        reference_layout=current_settings().reference_layout,
    )
    (out / "alleles.json").write_text(json.dumps(alleles_data, indent=2) + "\n")
    for key, fa in fastas.items():
        click.echo(f"{key}: {fa}")


@main.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Input consensus FASTA.",
)
@click.option(
    "--repeats-db",
    type=click.Path(exists=True),
    default=None,
    help="Custom repeat dictionary JSON.",
)
@click.option("--output-dir", "-o", type=click.Path(), default=".", help="Output directory.")
def classify(
    input_path: str,
    repeats_db: str | None,
    output_dir: str,
) -> None:
    """Classify repeat units in a consensus sequence."""
    from muc_one_span.classify import classify_sequence
    from muc_one_span.config import load_repeat_dictionary

    rd = load_repeat_dictionary(Path(repeats_db) if repeats_db else None)

    lines = Path(input_path).read_text().strip().splitlines()
    sequence = "".join(line for line in lines if not line.startswith(">"))

    result = classify_sequence(sequence, rd, settings=current_settings().classification)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    (out_dir / "repeats.json").write_text(json.dumps(result, indent=2) + "\n")
    (out_dir / "repeats.txt").write_text(result["structure"] + "\n")

    click.echo(f"Structure: {result['structure']}")
    if result["mutations_detected"]:
        click.echo(f"Mutations: {len(result['mutations_detected'])} detected")


@main.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(),
    help="Input FASTQ or BAM file.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(),
    default="results",
    help="Output directory.",
)
@click.option(
    "--reference",
    "-r",
    type=click.Path(),
    default=None,
    help="Reference FASTA (defaults to bundled ladder).",
)
@click.option(
    "--clair3-model",
    type=str,
    default=DEFAULT_SETTINGS.run.clair3_model,
    help="Path to Clair3 model.",
)
@click.option(
    "--threads", "-t", type=int, default=DEFAULT_SETTINGS.run.threads, help="Number of threads."
)
@click.option(
    "--min-coverage",
    type=int,
    default=DEFAULT_SETTINGS.run.min_coverage,
    help="Minimum read coverage.",
)
@click.option(
    "--min-qual",
    type=float,
    default=DEFAULT_SETTINGS.run.min_qual,
    help="Minimum QUAL score for VCF filtering (default 5.0).",
)
@click.option(
    "--report/--no-report",
    default=False,
    help="Generate HTML report (requires jinja2: pip install muc_one_span[report]).",
)
@click.option(
    "--report-igv",
    type=click.Choice(["embedded", "sidecar", "off"], case_sensitive=False),
    default="off",
    help="IGV alignment browser mode in HTML report (default: off).",
)
@click.option(
    "--platform",
    type=click.Choice(["hifi", "ont"], case_sensitive=False),
    default=DEFAULT_SETTINGS.run.platform,
    help="Sequencing platform (default: hifi).",
)
@click.option(
    "--mapping-timeout",
    type=float,
    default=DEFAULT_SETTINGS.run.mapping_timeout,
    help="Total mapping timeout in seconds (finite and positive; default 3600).",
)
@click.option(
    "--minimap2-preset",
    type=str,
    default=None,
    help="minimap2 -x preset (auto-selected from --platform if not set).",
)
@record_run_status
def run(
    input_path: str,
    output_dir: str,
    reference: str | None,
    clair3_model: str,
    threads: int,
    min_coverage: int,
    min_qual: float,
    report: bool,
    report_igv: str,
    platform: str,
    minimap2_preset: str | None,
    mapping_timeout: float = DEFAULT_SETTINGS.run.mapping_timeout,
) -> None:
    """Run the full MucOneSpan pipeline."""
    from muc_one_span.pipeline import execute_pipeline

    execute_pipeline(
        input_path,
        output_dir,
        reference,
        clair3_model,
        threads,
        min_coverage,
        min_qual,
        report,
        platform,
        minimap2_preset,
        report_igv=report_igv,
        mapping_timeout=mapping_timeout,
        settings=current_settings(),
        configuration=current_configuration_path(),
    )


@main.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to summary.json from a previous run.",
)
@click.option("--output", "-o", type=click.Path(), default="report.html", help="Output HTML path.")
@click.option("--sample-name", "-s", default=None, help="Sample name for report header.")
@click.option(
    "--repeats",
    type=click.Path(exists=True),
    default=None,
    help="Path to repeats.json for detailed repeat table.",
)
@click.option(
    "--report-igv",
    type=click.Choice(["embedded", "sidecar", "off"]),
    default="off",
    help="IGV report generation mode: embedded, sidecar, or off.",
)
@click.option(
    "--fasta", type=click.Path(exists=True), default=None, help="Reference FASTA for IGV."
)
@click.option("--bam", type=click.Path(exists=True), default=None, help="BAM file for IGV.")
@click.option("--vcf", type=click.Path(exists=True), default=None, help="VCF file for IGV.")
@click.option(
    "--allele-vcf",
    type=(str, click.Path(exists=True)),
    multiple=True,
    help="Labeled VCF track: LABEL PATH; repeat for each allele. Overrides --vcf.",
)
@click.option(
    "--run-status",
    type=click.Path(exists=True),
    default=None,
    help="Execution status JSON (default: summary sibling run_status.json).",
)
def report(
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
    """Generate an HTML report from pipeline results."""
    from muc_one_span.cli_report import execute_report

    execute_report(
        input_path,
        output,
        sample_name,
        repeats,
        report_igv,
        fasta,
        bam,
        vcf,
        allele_vcf,
        run_status,
    )


def _bundled_reference() -> Path:
    """Get path to the bundled reference ladder.

    Returns:
        Path to the bundled reference FASTA.

    Raises:
        SystemExit: If the bundled reference file does not exist.
    """
    import importlib.resources

    res = importlib.resources.files("muc_one_span.data.reference").joinpath("reference_ladder.fa")
    ref = Path(str(res))
    if not ref.exists():
        click.echo(
            f"Bundled reference not found at {ref}. Run 'muconespan ladder' to generate it.",
            err=True,
        )
        sys.exit(1)
    return ref


if __name__ == "__main__":
    main()
