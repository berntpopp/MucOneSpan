"""Read mapping with minimap2 and samtools."""

from __future__ import annotations

import logging
import subprocess as subprocess  # Retain the existing subprocess patch/import target.
import time
from pathlib import Path

from muc_one_span.settings import DEFAULT_SETTINGS
from muc_one_span.tool_pipeline import validate_timeout
from muc_one_span.tools import run_tool, run_tool_pipeline

PLATFORM_PRESETS: dict[str, str] = {"hifi": "map-hifi", "ont": "lr:hq"}
DEFAULT_MINIMAP2_PRESET = (
    DEFAULT_SETTINGS.run.minimap2_preset or PLATFORM_PRESETS[DEFAULT_SETTINGS.run.platform]
)

logger = logging.getLogger(__name__)


def bam_to_fastq(bam_path: Path, output_dir: Path) -> Path:
    """Convert BAM to FASTQ using samtools fastq.

    Captures stdout from ``samtools fastq`` and writes it to
    ``extracted_reads.fq`` in *output_dir*.

    Args:
        bam_path: Path to input BAM file.
        output_dir: Directory for output FASTQ.

    Returns:
        Path to the output FASTQ file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fastq_path = output_dir / "extracted_reads.fq"
    stdout = run_tool(["samtools", "fastq", str(bam_path)])
    fastq_path.write_text(stdout)
    return fastq_path


def map_reads(
    input_path: Path,
    reference_path: Path,
    output_dir: Path,
    threads: int = DEFAULT_SETTINGS.run.threads,
    preset: str = DEFAULT_MINIMAP2_PRESET,
    *,
    timeout: float = DEFAULT_SETTINGS.run.mapping_timeout,
) -> Path:
    """Map reads to reference using minimap2 and sort/index with samtools.

    Pipeline: ``minimap2 -a -x <preset>`` → ``samtools sort`` →
    ``samtools index``.  If *input_path* is a BAM file it is converted to
    FASTQ first with :func:`bam_to_fastq`.

    Args:
        input_path: Path to input FASTQ or BAM file.
        reference_path: Path to reference FASTA.
        output_dir: Directory for output files.
        threads: Number of threads for minimap2/samtools (default 4).
        preset: minimap2 preset passed via ``-x`` (default ``map-hifi``).
            Use ``lr:hq`` for Oxford Nanopore Q20+ reads.
        timeout: Finite positive seconds for conversion, mapping and indexing,
            including bounded process-group cleanup and diagnostic draining.

    Returns:
        Path to the sorted, indexed BAM file (``mapping.bam``).
    """
    validate_timeout(timeout)
    deadline = time.monotonic() + timeout
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Mapping reads from %s to %s", input_path.name, reference_path.name)

    bam_path = output_dir / "mapping.bam"
    artifacts = [
        bam_path,
        Path(str(bam_path) + ".bai"),
        Path(str(bam_path) + ".csi"),
        bam_path.with_suffix(".bai"),
        bam_path.with_suffix(".csi"),
    ]
    if input_path.resolve() in {artifact.resolve() for artifact in artifacts}:
        raise ValueError("Mapping input must differ from output BAM/index paths")
    try:
        for artifact in artifacts:
            artifact.unlink(missing_ok=True)
        actual_input = input_path
        if input_path.suffix.lower() == ".bam":
            actual_input = output_dir / "extracted_reads.fq"
            with actual_input.open("wb") as output:
                run_tool_pipeline(
                    [["samtools", "fastq", str(input_path)]],
                    timeout=_remaining(deadline),
                    stdout=output,
                )
        _run_mapping_pipeline(
            actual_input,
            reference_path,
            bam_path,
            threads,
            preset=preset,
            timeout=_remaining(deadline),
        )
        run_tool_pipeline([["samtools", "index", str(bam_path)]], timeout=_remaining(deadline))
    except BaseException:
        for artifact in artifacts:
            artifact.unlink(missing_ok=True)
        raise

    return bam_path


def _run_mapping_pipeline(
    input_path: Path,
    reference_path: Path,
    bam_path: Path,
    threads: int,
    preset: str = DEFAULT_MINIMAP2_PRESET,
    *,
    timeout: float = DEFAULT_SETTINGS.run.mapping_timeout,
) -> None:
    """Run minimap2 | samtools sort as a streaming pipeline.

    Pipes minimap2 SAM output directly into samtools sort, avoiding
    the need to hold the full SAM in memory.

    Args:
        input_path: Path to input FASTQ file.
        reference_path: Path to reference FASTA.
        bam_path: Path for sorted output BAM.
        threads: Number of threads for minimap2/samtools.
        preset: minimap2 preset passed via ``-x`` (default ``map-hifi``).

    Raises:
        FileNotFoundError: If minimap2 or samtools is not found.
        RuntimeError: If either process exits with non-zero status.
    """
    logger.info("Starting minimap2 | samtools sort pipeline (%d threads)", threads)
    minimap2_cmd = [
        "minimap2",
        "-a",
        "-x",
        preset,
        "-t",
        str(threads),
        str(reference_path),
        str(input_path),
    ]
    samtools_cmd = [
        "samtools",
        "sort",
        "-@",
        str(threads),
        "-o",
        str(bam_path),
    ]

    run_tool_pipeline([minimap2_cmd, samtools_cmd], timeout=timeout)


def _remaining(deadline: float) -> float:
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("Mapping timed out before the next stage")
    return remaining


def get_idxstats(bam_path: Path) -> str:
    """Run ``samtools idxstats`` and return the raw text output.

    Args:
        bam_path: Path to an indexed BAM file.

    Returns:
        Raw idxstats text output (tab-separated).
    """
    return run_tool(["samtools", "idxstats", str(bam_path)])
