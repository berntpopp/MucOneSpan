"""Display-only alignment of each allele's assigned reads for the hybrid IGV report.

Each allele shows exactly the reads the hybrid engine assigned to it (its spanning
members and the partial reads assigned by edit-distance competition): those reads
are written to a per-allele FASTQ and aligned to that allele's display contig only.
Aligning every read against both contigs instead lets near-identical alleles swap
reads, so the view would no longer match the evidence behind the call.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from muc_one_span.hybrid.igv_tracks import REFERENCE_PREFIX
from muc_one_span.hybrid.reads_io import read_input
from muc_one_span.mapping import map_reads
from muc_one_span.tools import run_tool


def _fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif line:
            records[name].append(line.strip())
    return {key: "".join(parts) for key, parts in records.items()}


def align_assigned_reads(
    input_path: Path,
    reference: Path,
    read_names: dict[str, list[str]],
    out_dir: Path,
    *,
    threads: int,
    preset: str,
    timeout: float,
) -> Path:
    """Align each allele's assigned reads to its own contig; return the merged BAM.

    ``read_names`` maps an allele (``allele_1``) to read names; its contig in
    ``reference`` is ``hybrid_<allele>``. Writes ``out_dir/<allele>/`` intermediates
    and ``out_dir/mapping.bam`` (sorted, indexed), indexes ``reference`` (``.fai``) for
    IGV Desktop, and removes the per-allele intermediates once merged.
    """
    contigs = _fasta_records(reference)
    owner = {name: allele for allele, names in read_names.items() for name in names}
    handles = {}
    try:
        for allele in read_names:
            (out_dir / allele).mkdir(parents=True, exist_ok=True)
            handles[allele] = (out_dir / allele / "reads.fastq").open("w", encoding="utf-8")
        for rec in read_input(input_path):
            if rec.name in owner:
                handles[owner[rec.name]].write(f"@{rec.name}\n{rec.seq}\n+\n{rec.qual}\n")
    finally:
        for handle in handles.values():
            handle.close()
    bams = []
    for allele in sorted(read_names):
        reads = out_dir / allele / "reads.fastq"
        if not read_names[allele] or reads.stat().st_size == 0:
            continue
        contig = f"{REFERENCE_PREFIX}{allele}"
        single = out_dir / allele / "reference.fa"
        single.write_text(f">{contig}\n{contigs[contig]}\n", encoding="utf-8")
        bams.append(
            map_reads(reads, single, out_dir / allele, threads, preset=preset, timeout=timeout)
        )
    if not bams:
        raise ValueError("IGV report: no assigned reads to display")
    merged = out_dir / "mapping.bam"
    run_tool(["samtools", "merge", "-f", "-o", str(merged), *map(str, bams)], timeout=timeout)
    run_tool(["samtools", "index", str(merged)], timeout=timeout)
    run_tool(["samtools", "faidx", str(reference)], timeout=timeout)
    for allele in read_names:
        shutil.rmtree(out_dir / allele, ignore_errors=True)
    return merged
