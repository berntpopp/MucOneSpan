"""Deterministic inputs for integration tests that invoke real genomics tools."""

from __future__ import annotations

import os
import random
import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def reference_data(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    rng = random.Random(1729)
    sequences = {name: "".join(rng.choices("ACGT", k=2500)) for name in ("contig_1", "contig_2")}
    reference = tmp_path / "reference.fa"
    reference.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
    return reference, sequences


@pytest.fixture
def aligned_bam(tmp_path: Path, reference_data: tuple[Path, dict[str, str]]) -> Path:
    """Three known primary alignments on two contigs, independent of minimap2."""
    _, sequences = reference_data
    sam = tmp_path / "input.sam"
    header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{len(seq)}\n" for name, seq in sequences.items()
    )
    records = []
    for name, contig, start in (
        ("read1", "contig_1", 100),
        ("read2", "contig_1", 500),
        ("read3", "contig_2", 200),
    ):
        sequence = sequences[contig][start : start + 1500]
        records.append(
            f"{name}\t0\t{contig}\t{start + 1}\t60\t1500M\t*\t0\t0\t{sequence}\t{'I' * 1500}\n"
        )
    sam.write_text(header + "".join(records))
    bam = tmp_path / "input.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    return bam


@pytest.fixture
def variant_data(
    tmp_path: Path, reference_data: tuple[Path, dict[str, str]]
) -> tuple[Path, dict[str, str]]:
    """SNP + insertion survive filtering; low QUAL and non-PASS records do not."""
    reference, sequences = reference_data
    subprocess.run(["samtools", "faidx", str(reference)], check=True)
    seq = sequences["contig_1"]
    alternate = next(base for base in "ACGT" if base != seq[499])
    low_alt = next(base for base in "ACGT" if base != seq[1499])
    fail_alt = next(base for base in "ACGT" if base != seq[1799])
    vcf = tmp_path / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        + "".join(f"##contig=<ID={name},length={len(s)}>\n" for name, s in sequences.items())
        + '##FILTER=<ID=LowQual,Description="Failed caller filter">\n'
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        + f"contig_1\t500\t.\t{seq[499]}\t{alternate}\t30\tPASS\t.\tGT\t1/1\n"
        + f"contig_1\t1000\t.\t{seq[999]}\t{seq[999]}AC\t30\tPASS\t.\tGT\t1/1\n"
        + f"contig_1\t1500\t.\t{seq[1499]}\t{low_alt}\t2\tPASS\t.\tGT\t1/1\n"
        + f"contig_1\t1800\t.\t{seq[1799]}\t{fail_alt}\t30\tLowQual\t.\tGT\t1/1\n"
    )
    expected = {**sequences, "contig_1": seq[:499] + alternate + seq[500:1000] + "AC" + seq[1000:]}
    return vcf, expected


@pytest.fixture
def clair3_model() -> Path:
    """Use an explicit model or Clair3's adjacent bundled HiFi model."""
    configured = os.environ.get("CLAIR3_MODEL")
    executable = shutil.which("run_clair3.sh")
    if configured:
        model = Path(configured)
        if not model.is_dir():
            pytest.fail(f"CLAIR3_MODEL directory does not exist: {model}")
    elif executable:
        model = Path(executable).resolve().parent / "models" / "hifi"
    else:
        pytest.skip("Clair3 not installed (run_clair3.sh is absent from PATH)")
    if not all((model / f"{name}.index").is_file() for name in ("pileup", "full_alignment")):
        pytest.skip("Clair3 HiFi model checkpoints unavailable; set CLAIR3_MODEL")
    return model
