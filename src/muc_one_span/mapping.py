"""Read mapping with minimap2 and samtools."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from muc_one_span.settings import DEFAULT_SETTINGS
from muc_one_span.tools import run_tool

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

    Returns:
        Path to the sorted, indexed BAM file (``mapping.bam``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Mapping reads from %s to %s", input_path.name, reference_path.name)

    # If input is BAM, convert to FASTQ first
    actual_input = input_path
    if input_path.suffix.lower() == ".bam":
        actual_input = bam_to_fastq(input_path, output_dir)

    bam_path = output_dir / "mapping.bam"

    # Pipe minimap2 SAM output directly to samtools sort to avoid
    # holding the full SAM in memory. The pipeline implementation uses
    # subprocesses with argument lists (no shell invocation).
    _run_mapping_pipeline(actual_input, reference_path, bam_path, threads, preset=preset)

    # samtools index
    run_tool(["samtools", "index", str(bam_path)])

    return bam_path


def _run_mapping_pipeline(
    input_path: Path,
    reference_path: Path,
    bam_path: Path,
    threads: int,
    preset: str = DEFAULT_MINIMAP2_PRESET,
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

    try:
        p1 = subprocess.Popen(
            minimap2_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError("Tool not found: minimap2") from exc

    # Verify stdout was captured before piping into samtools
    if p1.stdout is None:
        p1.kill()
        p1.wait()
        raise RuntimeError("minimap2 process stdout was not captured")

    try:
        p2 = subprocess.Popen(
            samtools_cmd,
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError as exc:
        p1.kill()
        p1.wait()
        raise FileNotFoundError("Tool not found: samtools") from exc

    # Allow p1 to receive SIGPIPE if p2 exits early
    p1.stdout.close()

    # Drain p1.stderr concurrently with p2.communicate() to avoid deadlock.
    # If we read p1.stderr first, p2 could fill its stderr buffer and block;
    # if we call p2.communicate() first, p1 could fill its stderr buffer.
    import threading

    p1_stderr_chunks: list[bytes] = []

    def _drain_p1_stderr() -> None:
        if p1.stderr:
            p1_stderr_chunks.append(p1.stderr.read())

    stderr_thread = threading.Thread(target=_drain_p1_stderr)
    stderr_thread.start()
    _, p2_stderr = p2.communicate()
    stderr_thread.join()
    p1.wait()

    p1_stderr = p1_stderr_chunks[0].decode() if p1_stderr_chunks else ""

    if p1.returncode != 0:
        raise RuntimeError(f"minimap2 failed with exit code {p1.returncode}.\nstderr: {p1_stderr}")
    if p2.returncode != 0:
        raise RuntimeError(
            f"samtools sort failed with exit code {p2.returncode}.\nstderr: {p2_stderr.decode()}"
        )


def get_idxstats(bam_path: Path) -> str:
    """Run ``samtools idxstats`` and return the raw text output.

    Args:
        bam_path: Path to an indexed BAM file.

    Returns:
        Raw idxstats text output (tab-separated).
    """
    return run_tool(["samtools", "idxstats", str(bam_path)])
