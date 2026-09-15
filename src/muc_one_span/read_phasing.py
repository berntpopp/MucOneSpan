"""Optional read-backed phase for unresolved same-length candidates.

Phase links retained variants; it does not verify reference-filled bases, allele
lengths, or biological molecule independence. WhatsHap is optional and its
installed version may leave multiallelic sites unresolved.
"""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import tempfile
from pathlib import Path

from muc_one_span.phasing import phase_evidence
from muc_one_span.settings import DEFAULT_SETTINGS, ReadPhasingSettings
from muc_one_span.tools import _clean_path_for_externals, run_tool, run_tool_iter
from muc_one_span.vcf import parse_vcf_genotypes, select_vcf_sample

logger = logging.getLogger(__name__)


def _primary_records(bam_path: Path, sam_path: Path, source_map: Path) -> dict:
    """Stream primary records, changing QNAME only and preserving source order.

    Input ordinals are one-based alignment-record positions, including excluded
    records and excluding headers. Different primary records remain independent
    observations even when their original names collide.
    """
    total, retained = 0, 0
    with sam_path.open("w") as sam, source_map.open("w") as mapping:
        for line in run_tool_iter(["samtools", "view", "-h", str(bam_path)]):
            if line.startswith("@"):
                sam.write(line)
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("Malformed alignment record in read-phasing input")
            flag = int(fields[1])
            if flag & (4 | 256 | 2048):
                continue
            phase_name = f"phase_record_{retained:012d}"
            mapping.write(
                json.dumps(
                    {
                        "input_record_ordinal": total,
                        "original_qname": fields[0],
                        "phase_qname": phase_name,
                    }
                )
                + "\n"
            )
            fields[0] = phase_name
            sam.write("\t".join(fields) + "\n")
            retained += 1
    return {
        "input_alignment_records": total,
        "primary_records": retained,
        "excluded_records": total - retained,
    }


def _variant_identity(variants: list[dict]) -> list[tuple]:
    """Compare selected biological alleles, allowing whole-haplotype swaps."""
    return [
        (
            variant["chrom"],
            variant["pos"],
            variant["ref"],
            variant["alt"],
            tuple(sorted(re.split(r"[/|]", variant["genotype"]))),
        )
        for variant in variants
    ]


def _assigned_read_counts(
    read_list: Path,
    primary_records: int,
    sample: str,
    variants: list[dict],
    phased: bool,
) -> tuple[int, dict[str, dict[str, int]]]:
    """Validate tool assignments against the selected sample, block and records.

    Counts describe the internal selection only, grouped by phase set and the
    tool's zero-based haplotype labels. A nonempty coherent read list is required
    for a phase claim; it is not a minimum biological support rule.
    """
    heterozygous = [
        variant for variant in variants if len(set(re.split(r"[/|]", variant["genotype"]))) == 2
    ]
    positions = {variant["pos"] for variant in heterozygous}
    phase_sets = {str(variant.get("phase_set")) for variant in heterozygous}
    header = (
        "#readname\tsource_id\tsample\tphaseset\thaplotype\tcovered_variants"
        "\tfirst_variant_pos\tlast_variant_pos"
    )
    assigned: set[str] = set()
    selected_counts: dict[str, dict[str, int]] = {}
    with read_list.open() as handle:
        if handle.readline().rstrip("\n") != header:
            raise ValueError("Malformed phasing read-list header")
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                raise ValueError("Malformed phasing read assignment")
            name, source, selected, phase_set, haplotype, count, first, last = fields
            match = re.fullmatch(r"phase_record_([0-9]{12,})", name)
            if match is None or int(match[1]) >= primary_records or name in assigned:
                raise ValueError("Invalid or duplicate read identity in phasing read list")
            if (
                source != "0"
                or selected != sample
                or haplotype not in {"0", "1"}
                or not all(value.isdecimal() for value in (phase_set, count, first, last))
            ):
                raise ValueError("Invalid sample, source or values in phasing read list")
            covered, start, end = int(count), int(first), int(last)
            if (
                int(phase_set) < 1
                or start not in positions
                or end not in positions
                or not 0 < covered <= sum(start <= pos <= end for pos in positions)
                or (phased and phase_set not in phase_sets)
            ):
                raise ValueError("Contradictory phase set or variant count in phasing read list")
            assigned.add(name)
            selected_counts.setdefault(phase_set, {"0": 0, "1": 0})[haplotype] += 1
    if phased and not assigned:
        raise ValueError("Phased output has no phasing read assignments")
    return len(assigned), selected_counts


def phase_same_length_reads(
    vcf_path: Path,
    bam_path: Path,
    reference_path: Path,
    output_dir: Path,
    *,
    sample: str | None = None,
    settings: ReadPhasingSettings | None = None,
) -> tuple[Path, dict]:
    """Optionally resolve multiple heterozygous loci with existing remapped reads.

    Return a selected VCF and additive ``read_phasing`` metadata. Already usable,
    empty, incomplete, non-diploid and conflicting evidence is preserved. Missing
    WhatsHap or absent primary reads returns explicit unresolved metadata; a
    present tool's execution failure propagates. Malformed VCF input also fails.

    Only a result accepted by the existing common-phase-set policy replaces the
    input. Other attempts retain the original VCF and expose their candidate
    status. Every attempt gets its own directory to avoid stale-output reuse.
    Temporary alignment copies are removed; the source map and tool read list
    persist. Read-list counts reflect internal selection/downsampling, not all BAM
    records. Per-phase-set counts retain the tool's zero-based haplotype labels;
    zero selected records is not evidence that a biological haplotype is absent.
    """
    settings = settings or DEFAULT_SETTINGS.read_phasing
    selected_sample = select_vcf_sample(vcf_path, sample)
    variants = parse_vcf_genotypes(vcf_path, sample=selected_sample)
    original_status = phase_evidence(variants)["phase_status"]
    metadata: dict = {
        "method": "whatshap",
        "version": None,
        "status": "not_needed",
        "sample": selected_sample,
        "input_phase_status": original_status,
        "output_phase_status": original_status,
        "primary_records": None,
        "source_map": None,
        "read_list": None,
        "phase_command": None,
        "tool_parameter_overrides": {
            "internal_downsampling": settings.internal_downsampling,
            "mapping_quality": settings.mapping_quality,
        },
        "unspecified_parameters": "installed_tool_defaults; version recorded when invoked",
        "phasing_selected_records_by_phase_set": None,
        "phasing_read_count_scope": "internal_selection_not_all_primary_records",
    }
    if original_status not in {"unphased", "missing_phase_set", "disconnected_phase_sets"}:
        return vcf_path, metadata
    external_path = _clean_path_for_externals(os.environ.get("PATH", ""))
    if shutil.which("whatshap", path=external_path) is None:
        metadata["status"] = "unavailable"
        return vcf_path, metadata

    metadata["version"] = run_tool(["whatshap", "--version"]).strip()
    output_dir.mkdir(parents=True, exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix="attempt_", dir=output_dir))
    source_map = attempt / "read_sources.jsonl"
    read_list = attempt / "assigned_reads.tsv"
    candidate = attempt / "candidate.vcf"
    metadata.update(source_map=str(source_map), read_list=str(read_list))
    with tempfile.TemporaryDirectory(prefix="alignments_", dir=attempt) as temporary:
        sam_path = Path(temporary) / "primary.sam"
        metadata.update(_primary_records(bam_path, sam_path, source_map))
        if not metadata["primary_records"]:
            metadata.update(status="insufficient_reads", read_list=None)
            return vcf_path, metadata
        phase_bam = Path(temporary) / "primary.bam"
        run_tool(["samtools", "view", "-b", "-o", str(phase_bam), str(sam_path)])
        run_tool(["samtools", "index", str(phase_bam)])
        overrides: list[str] = []
        if settings.internal_downsampling is not None:
            overrides.extend(["--internal-downsampling", str(settings.internal_downsampling)])
        if settings.mapping_quality is not None:
            overrides.extend(["--mapping-quality", str(settings.mapping_quality)])
        phase_command = [
            "whatshap",
            "phase",
            "--indels",
            "--ignore-read-groups",
            *overrides,
            "--sample",
            selected_sample,
            "--reference",
            str(reference_path),
            "--output-read-list",
            str(read_list),
            "-o",
            str(candidate),
            str(vcf_path),
            str(phase_bam),
        ]
        metadata["phase_command"] = phase_command
        run_tool(phase_command)
    phased_variants = parse_vcf_genotypes(candidate, sample=selected_sample)
    if _variant_identity(variants) != _variant_identity(phased_variants):
        raise ValueError("Read phaser changed variant identity or genotype")
    candidate_status = phase_evidence(phased_variants)["phase_status"]
    assigned_records, selected_counts = _assigned_read_counts(
        read_list,
        metadata["primary_records"],
        selected_sample,
        phased_variants,
        candidate_status == "phased",
    )
    metadata.update(
        phasing_assigned_records=assigned_records,
        phasing_selected_records_by_phase_set=selected_counts,
        candidate_vcf=str(candidate),
        candidate_phase_status=candidate_status,
    )
    if candidate_status != "phased":
        metadata["status"] = "unresolved"
        return vcf_path, metadata

    phased_vcf = attempt / "phased.vcf.gz"
    run_tool(["bcftools", "view", "-Oz", "-o", str(phased_vcf), str(candidate)])
    run_tool(["bcftools", "index", str(phased_vcf)])
    metadata.update(status="phased", output_phase_status="phased", selected_vcf=str(phased_vcf))
    return phased_vcf, metadata


def haplotag_and_split_reads(
    phased_vcf: Path,
    bam_path: Path,
    reference_path: Path,
    output_dir: Path,
    *,
    min_reads: int = 5,
) -> tuple[Path, Path, int, int] | None:
    """Haplotag an alignment BAM using phased variants and split into two haplotype BAMs.

    Args:
        phased_vcf: VCF with phased heterozygous variants and phase set.
        bam_path: Input BAM file (e.g. remapped cluster reads).
        reference_path: Contig reference FASTA.
        output_dir: Destination directory.
        min_reads: Minimum reads required per haplotype to proceed with separation.

    Returns:
        (hp1_bam, hp2_bam, hp1_count, hp2_count) if both haplotypes have >= min_reads,
        or None if haplotagging was unresolvable or had insufficient reads.
    """
    external_path = _clean_path_for_externals(os.environ.get("PATH", ""))
    if shutil.which("whatshap", path=external_path) is None:
        logger.info("WhatsHap unavailable for haplotagging; skipping BAM split")
        return None

    output_dir.mkdir(parents=True, exist_ok=True)
    haplotagged_bam = output_dir / "haplotagged.bam"
    try:
        run_tool(
            [
                "whatshap",
                "haplotag",
                "--ignore-read-groups",
                "--reference",
                str(reference_path),
                "-o",
                str(haplotagged_bam),
                str(phased_vcf),
                str(bam_path),
            ]
        )
    except Exception as exc:
        logger.warning("WhatsHap haplotagging failed: %s", exc)
        return None

    hp1_bam = output_dir / "allele_1.bam"
    hp2_bam = output_dir / "allele_2.bam"
    try:
        run_tool(["samtools", "view", "-b", "-d", "HP:1", "-o", str(hp1_bam), str(haplotagged_bam)])
        run_tool(["samtools", "index", str(hp1_bam)])
        run_tool(["samtools", "view", "-b", "-d", "HP:2", "-o", str(hp2_bam), str(haplotagged_bam)])
        run_tool(["samtools", "index", str(hp2_bam)])

        c1_str = run_tool(["samtools", "view", "-c", str(hp1_bam)]).strip()
        c2_str = run_tool(["samtools", "view", "-c", str(hp2_bam)]).strip()
        c1 = int(c1_str) if c1_str.isdigit() else 0
        c2 = int(c2_str) if c2_str.isdigit() else 0
    except Exception as exc:
        logger.warning("Failed to extract haplotype BAMs: %s", exc)
        return None
    finally:
        haplotagged_bam.unlink(missing_ok=True)

    if c1 < min_reads or c2 < min_reads:
        logger.info(
            "Insufficient reads in haplotagged partitions: hp1=%d, hp2=%d (min %d)",
            c1,
            c2,
            min_reads,
        )
        return None

    logger.info("Successfully haplotagged and split reads: hp1=%d, hp2=%d", c1, c2)
    return hp1_bam, hp2_bam, c1, c2
