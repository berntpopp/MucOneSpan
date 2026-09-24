"""The ``run`` command: full pipeline execution (moved from ``cli.py``)."""

from __future__ import annotations

import click

from muc_one_span.cli_settings import current_configuration_path, current_settings
from muc_one_span.run_status import record_run_status
from muc_one_span.settings import DEFAULT_SETTINGS


@click.command()
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
    help="Minimum VCF QUAL (default 5.0); see calling.haploid_min_qual for length-split calls.",
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
@click.option(
    "--engine",
    type=click.Choice(["ladder", "hybrid"]),
    default=DEFAULT_SETTINGS.run.engine,
    help="Allele reconstruction engine (hybrid is experimental; default: ladder).",
)
@click.option(
    "--assay",
    type=click.Choice(["amplicon", "genomic"]),
    default=DEFAULT_SETTINGS.run.assay,
    help="Library type used by the hybrid engine (default: amplicon).",
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
    engine: str = DEFAULT_SETTINGS.run.engine,
    assay: str = DEFAULT_SETTINGS.run.assay,
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
        engine=engine,
        assay=assay,
    )
