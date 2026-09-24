"""Variant calling with Clair3 and VCF processing with bcftools."""

from __future__ import annotations

import copy
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from muc_one_span.mapping import DEFAULT_MINIMAP2_PRESET
from muc_one_span.phasing import (
    annotate_consensus_candidate,
    length_partition_selection,
    phase_evidence,
)
from muc_one_span.read_phasing import haplotag_and_split_reads, phase_same_length_reads
from muc_one_span.settings import DEFAULT_SETTINGS, CallingSettings, ReadPhasingSettings
from muc_one_span.tools import run_tool
from muc_one_span.vcf import filter_vcf, parse_vcf_genotypes, parse_vcf_variants

__all__ = [
    "call_variants_per_allele",
    "disambiguate_same_length_alleles",
    "extract_allele_reads",
    "filter_vcf",
    "parse_vcf_genotypes",
    "parse_vcf_variants",
    "run_clair3",
]

logger = logging.getLogger(__name__)


def _haploid_filter(
    min_qual: float, settings: CallingSettings
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return length-partitioned VCF filter options and their explicit provenance.

    ``calling.haploid_min_qual`` (default 4.0) applies to haploid, length-partitioned
    calls; ``null`` makes them follow ``run.min_qual``/``--min-qual``. An explicit
    non-default ``--min-qual`` that is overridden is logged, never silently dropped.
    """
    if settings.haploid_min_qual is None:
        qual, source = float(min_qual), "run.min_qual"
    else:
        qual, source = float(settings.haploid_min_qual), "calling.haploid_min_qual"
        if min_qual != DEFAULT_SETTINGS.run.min_qual:
            logger.warning(
                "--min-qual %s is not applied to length-partitioned calls; "
                "calling.haploid_min_qual=%s is used. Set it to null to follow --min-qual.",
                min_qual,
                qual,
            )
    options = {
        "haploid_majority": settings.haploid_majority,
        "haploid_min_qual": qual,
        "haploid_alt_fraction": settings.haploid_alt_fraction,
        "haploid_ref_fraction": settings.haploid_ref_fraction,
    }
    provenance = {
        "min_qual": qual,
        "min_qual_source": source,
        "haploid_majority": settings.haploid_majority,
    }
    return options, provenance


def extract_allele_reads(
    bam_path: Path,
    contig_names: str | list[str],
    output_dir: Path,
    *,
    output_name: str = "allele_reads.bam",
) -> Path:
    """Extract reads mapped to one or more contigs from a BAM file.

    Uses ``samtools view -b`` to subset the BAM to the specified contig(s)
    and indexes the result.

    Args:
        bam_path: Path to the full mapping BAM.
        contig_names: Single contig name or list of contig names to extract.
        output_dir: Directory for output files.
        output_name: File name of the extracted BAM inside ``output_dir``.

    Returns:
        Path to the extracted, indexed BAM file.
    """
    if isinstance(contig_names, str):
        contig_names = [contig_names]

    output_dir.mkdir(parents=True, exist_ok=True)
    out_bam = output_dir / output_name

    run_tool(
        [
            "samtools",
            "view",
            "-b",
            "-o",
            str(out_bam),
            str(bam_path),
            *contig_names,
        ]
    )
    run_tool(["samtools", "index", str(out_bam)])

    return out_bam


def _extract_and_remap_reads(
    bam_path: Path,
    cluster_contigs: list[str],
    peak_contig: str,
    reference_path: Path,
    output_dir: Path,
    threads: int = DEFAULT_SETTINGS.run.threads,
    preset: str = DEFAULT_MINIMAP2_PRESET,
) -> Path:
    """Extract reads from cluster contigs, convert to FASTQ, remap to peak contig.

    Reads from the ladder mapping are spread across multiple contigs in a
    cluster (e.g. contig_48 through contig_54).  Clair3 needs all reads
    aligned to a *single* reference contig to call variants.  This function
    extracts the cluster reads, converts to FASTQ, extracts the single peak
    contig as a mini-reference, and remaps with minimap2.

    Args:
        bam_path: Full ladder mapping BAM.
        cluster_contigs: All contig names in the allele's cluster.
        peak_contig: The single peak contig name to remap against.
        reference_path: Full ladder reference FASTA (for extracting the contig).
        output_dir: Working directory for intermediate files.
        threads: Thread count for minimap2/samtools.
        preset: minimap2 preset (default ``"map-hifi"``). Use ``"lr:hq"`` for ONT.

    Returns:
        Path to the remapped, sorted, indexed BAM.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Extract reads from all cluster contigs (ladder coordinates; intermediate)
    cluster_bam = extract_allele_reads(
        bam_path, cluster_contigs, output_dir, output_name="cluster_reads.bam"
    )

    # 2. Convert to FASTQ, then drop the cluster BAM so allele_reads.bam is only
    #    ever the remapped single-contig BAM.
    fastq_path = output_dir / "cluster_reads.fq"
    stdout = run_tool(["samtools", "fastq", str(cluster_bam)])
    fastq_path.write_text(stdout)
    cluster_bam.unlink(missing_ok=True)
    Path(f"{cluster_bam}.bai").unlink(missing_ok=True)

    # 3. Extract peak contig as mini-reference
    contig_ref = output_dir / f"{peak_contig}.fa"
    stdout = run_tool(
        [
            "samtools",
            "faidx",
            str(reference_path),
            peak_contig,
        ]
    )
    contig_ref.write_text(stdout)
    run_tool(["samtools", "faidx", str(contig_ref)])

    # 4. Remap to peak contig
    sam_path = output_dir / "remapped.sam"
    sam_output = run_tool(
        [
            "minimap2",
            "-a",
            "-x",
            preset,
            "-t",
            str(threads),
            str(contig_ref),
            str(fastq_path),
        ]
    )
    sam_path.write_text(sam_output)

    # 5. Sort and index
    remapped_bam = output_dir / "allele_reads.bam"
    run_tool(
        [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(remapped_bam),
            str(sam_path),
        ]
    )
    run_tool(["samtools", "index", str(remapped_bam)])

    # Clean up intermediates
    sam_path.unlink(missing_ok=True)
    fastq_path.unlink(missing_ok=True)

    return remapped_bam


def run_clair3(
    bam_path: Path,
    reference_path: Path,
    output_dir: Path,
    model_path: str = DEFAULT_SETTINGS.run.clair3_model,
    platform: str = DEFAULT_SETTINGS.run.platform,
    threads: int = DEFAULT_SETTINGS.run.threads,
    *,
    settings: CallingSettings | None = None,
) -> Path:
    """Run the Clair3 variant caller on a BAM file.

    Args:
        bam_path: Path to input BAM (reads for one allele).
        reference_path: Path to the reference FASTA.
        output_dir: Directory for Clair3 output.
        model_path: Path to Clair3 model directory (optional).
        platform: Sequencing platform (default ``"hifi"``).
        threads: Number of threads (default 4).

    Returns:
        Path to the Clair3 output VCF (``merge_output.vcf.gz``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Running Clair3 on %s", bam_path.name)
    settings = settings or DEFAULT_SETTINGS.calling

    cmd = [
        "run_clair3.sh",
        f"--bam_fn={bam_path}",
        f"--ref_fn={reference_path}",
        f"--output={output_dir}",
        f"--threads={threads}",
        f"--platform={platform}",
        f"--sample_name={settings.sample_name}",
        # Our contigs are named contig_N, not chr1..22/X/Y, so Clair3
        # must be told to process all contigs.
        "--include_all_ctgs",
    ]
    if model_path:
        cmd.append(f"--model_path={model_path}")

    run_tool(cmd)

    return output_dir / "merge_output.vcf.gz"


def disambiguate_same_length_alleles(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = DEFAULT_SETTINGS.run.clair3_model,
    threads: int = DEFAULT_SETTINGS.run.threads,
    min_qual: float = float(DEFAULT_SETTINGS.run.min_qual),
    min_dp: int = 5,
    platform: str = DEFAULT_SETTINGS.run.platform,
    preset: str = DEFAULT_MINIMAP2_PRESET,
    *,
    read_phase: bool | None = None,
    settings: CallingSettings | None = None,
    read_phasing_settings: ReadPhasingSettings | None = None,
) -> dict[str, Path]:
    """Disambiguate same-length alleles using Clair3 genotype calls.

    ``read_phase=True`` enables experimental read-backed phase, disabled by
    default after an extra mutation call in cached development validation.

    Runs Clair3 on the cluster reads and interprets complete genotype indices.
    Cross-site phase requires a common phase set; a single heterozygous site
    provides an unordered pair. Missing or disconnected evidence stays mixed.

    Returns only allele-key to VCF-path mappings. Allele dictionaries receive
    explicit phase and sequence-evidence status. Unresolved multiple sites keep
    one mixed candidate; supported haplotypes share a VCF and select GT indices.
    """
    settings = settings or DEFAULT_SETTINGS.calling
    read_phase = settings.read_phase if read_phase is None else read_phase
    allele_info = alleles["allele_1"]
    contig_name = allele_info["contig_name"]
    cluster_contigs = allele_info["cluster_contigs"]
    merged_dir = output_dir / "merged"

    # Remap all cluster reads to peak contig
    merged_bam = _extract_and_remap_reads(
        bam_path,
        cluster_contigs,
        contig_name,
        reference_path,
        merged_dir,
        threads,
        preset=preset,
    )

    # Run Clair3
    contig_ref = merged_dir / f"{contig_name}.fa"
    clair3_dir = merged_dir / "clair3"
    raw_vcf = run_clair3(
        merged_bam,
        contig_ref,
        clair3_dir,
        model_path=clair3_model,
        platform=platform,
        threads=threads,
        settings=settings,
    )
    filtered_vcf = filter_vcf(raw_vcf, contig_ref, merged_dir, min_qual=min_qual, min_dp=min_dp)

    read_phasing = {"status": "experimental_disabled", "method": "whatshap"}
    if read_phase:
        filtered_vcf, read_phasing = phase_same_length_reads(
            filtered_vcf,
            merged_bam,
            contig_ref,
            merged_dir / "read_phasing",
            settings=read_phasing_settings,
        )

    if read_phasing.get("output_phase_status") == "phased":
        split_result = haplotag_and_split_reads(
            filtered_vcf,
            merged_bam,
            contig_ref,
            merged_dir / "haplotag",
            min_reads=min_dp,
        )
        if split_result is not None:
            hp1_bam, hp2_bam, hp1_count, hp2_count = split_result
            per_allele_threads = max(1, threads // 2)
            hp_filter, _ = _haploid_filter(min_qual, settings)

            def _call_hp(hp_key: str, hp_bam: Path) -> tuple[str, Path]:
                hp_dir = merged_dir / hp_key
                vcf_raw = run_clair3(
                    hp_bam,
                    contig_ref,
                    hp_dir / "clair3",
                    model_path=clair3_model,
                    platform=platform,
                    threads=per_allele_threads,
                    settings=settings,
                )
                vcf_filtered = filter_vcf(
                    vcf_raw,
                    contig_ref,
                    hp_dir,
                    min_qual=min_qual,
                    min_dp=min_dp,
                    **hp_filter,
                )
                return hp_key, vcf_filtered

            with ThreadPoolExecutor(max_workers=2) as executor:
                hp_futures = [
                    executor.submit(_call_hp, "allele_1", hp1_bam),
                    executor.submit(_call_hp, "allele_2", hp2_bam),
                ]
                hp_results = dict(f.result() for f in hp_futures)

            v_1 = parse_vcf_genotypes(hp_results["allele_1"])
            v_2 = parse_vcf_genotypes(hp_results["allele_2"])

            calls_1 = {(v["chrom"], v["pos"], v["ref"], v["alt"]): v["genotype"] for v in v_1}
            calls_2 = {(v["chrom"], v["pos"], v["ref"], v["alt"]): v["genotype"] for v in v_2}
            is_distinct = calls_1 != calls_2

            alleles["homozygous"] = not is_distinct
            alleles["sequence_identity_status"] = (
                "resolved_distinct" if is_distinct else "unresolved"
            )
            alleles["phase_status"] = "phased"

            ev_1 = phase_evidence(v_1)
            s_1 = v_1[0].get("sample") if v_1 else None
            annotate_consensus_candidate(
                alleles["allele_1"], ev_1, 1, s_1, str(hp_results["allele_1"])
            )
            alleles["allele_1"]["read_phasing"] = read_phasing
            alleles["allele_1"]["reads"] = hp1_count
            alleles["allele_1"]["independent_haplotype_evidence"] = True
            alleles["allele_1"]["sequence_identity_status"] = (
                "resolved_distinct" if is_distinct else "unresolved"
            )

            if "allele_2" not in alleles:
                alleles["allele_2"] = copy.deepcopy(allele_info)
            alleles["allele_2"].pop("candidate_duplicate_of", None)
            ev_2 = phase_evidence(v_2)
            s_2 = v_2[0].get("sample") if v_2 else None
            annotate_consensus_candidate(
                alleles["allele_2"], ev_2, 1, s_2, str(hp_results["allele_2"])
            )
            alleles["allele_2"]["read_phasing"] = read_phasing
            alleles["allele_2"]["reads"] = hp2_count
            alleles["allele_2"]["independent_haplotype_evidence"] = True
            alleles["allele_2"]["sequence_identity_status"] = (
                "resolved_distinct" if is_distinct else "unresolved"
            )

            return hp_results

    variants = parse_vcf_genotypes(filtered_vcf)
    evidence = phase_evidence(variants)
    # A single observed length is not evidence of sequence identity.
    alleles["homozygous"] = False
    alleles["sequence_identity_status"] = "unresolved"
    alleles["phase_status"] = evidence["phase_status"]
    sample = variants[0].get("sample") if variants else None
    results: dict[str, Path] = {}
    for index, haplotype in enumerate(evidence["haplotypes"], start=1):
        key = f"allele_{index}"
        if key not in alleles:
            alleles[key] = copy.deepcopy(allele_info)
        if haplotype in (1, 2):
            alleles[key].pop("candidate_duplicate_of", None)
        annotate_consensus_candidate(alleles[key], evidence, haplotype, sample, str(filtered_vcf))
        alleles[key]["read_phasing"] = read_phasing
        results[key] = filtered_vcf
    if "allele_2" in alleles and "allele_2" not in results:
        alias = alleles["allele_2"]
        alias.update({key: value for key, value in evidence.items() if key != "haplotypes"})
        alias.update(
            reconstruction_status="not_separately_resolved",
            candidate_duplicate_of="allele_1",
            independent_haplotype_evidence=False,
        )
        for field in (
            "vcf_path",
            "consensus_haplotype",
            "consensus_sample",
            "sequence_source",
            "variant_observation_group",
            "consensus_context",
            "consensus_policy",
            "read_phasing",
            "vntr_phase_status",
        ):
            alias.pop(field, None)
    return results


def call_variants_per_allele(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = DEFAULT_SETTINGS.run.clair3_model,
    threads: int = DEFAULT_SETTINGS.run.threads,
    min_qual: float = float(DEFAULT_SETTINGS.run.min_qual),
    min_dp: int = 5,
    platform: str = DEFAULT_SETTINGS.run.platform,
    preset: str = DEFAULT_MINIMAP2_PRESET,
    *,
    read_phase: bool | None = None,
    settings: CallingSettings | None = None,
    read_phasing_settings: ReadPhasingSettings | None = None,
) -> dict[str, Path]:
    """Run variant calling for each detected allele.

    For each allele key in *alleles* (``"allele_1"`` and optionally
    ``"allele_2"``), this function:

    1. Extracts the reads mapped to the allele's best-matching contig.
    2. Calls variants with Clair3.
    3. Filters the VCF with :func:`filter_vcf`.

    For same-length alleles, delegates to :func:`disambiguate_same_length_alleles`.

    Args:
        bam_path: Path to the full mapping BAM.
        reference_path: Path to the ladder reference FASTA.
        alleles: Allele detection result from :func:`~muc_one_span.alleles.detect_alleles`.
        output_dir: Base output directory.
        clair3_model: Path to Clair3 model directory (optional).
        threads: Number of threads (default 4).
        min_qual: Minimum QUAL score for VCF filtering (default 5.0).
            Low-confidence variants are retained but penalized in the
            confidence scoring system rather than hard-filtered.
        min_dp: Accepted for API compatibility; not applied by VCF filtering.
        platform: Sequencing platform for Clair3 (default ``"hifi"``).
        preset: minimap2 preset for remapping (default ``"map-hifi"``).
        read_phase: Experimental library-only read phasing, disabled by default.

    Returns:
        Dictionary mapping allele key (``"allele_1"`` / ``"allele_2"``) to
        the filtered VCF path.
    """
    settings = settings or DEFAULT_SETTINGS.calling
    read_phase = settings.read_phase if read_phase is None else read_phase
    if alleles.get("same_length"):
        logger.info("Same-length alleles detected, using disambiguation")
        disambig = disambiguate_same_length_alleles(
            bam_path,
            reference_path,
            alleles,
            output_dir,
            clair3_model,
            threads,
            min_qual,
            min_dp,
            platform=platform,
            preset=preset,
            read_phase=read_phase,
            settings=settings,
            read_phasing_settings=read_phasing_settings,
        )
        return disambig

    # Collect allele keys to process (skip allele_2 when homozygous)
    allele_keys = [
        k
        for k in ("allele_1", "allele_2")
        if k in alleles and not (alleles.get("homozygous") and k == "allele_2")
    ]

    filter_options, filter_provenance = _haploid_filter(min_qual, settings)

    def _process_allele(allele_key: str) -> tuple[str, Path]:
        allele_info = alleles[allele_key]
        # Use the peak contig name from allele detection (contig_N where N is
        # canonical X repeat count, NOT total length).  Fall back to computing
        # from length for backwards compatibility with older alleles.json.
        contig_name = allele_info.get("contig_name", f"contig_{allele_info['length']}")
        cluster_contigs = allele_info.get("cluster_contigs", [contig_name])
        allele_dir = output_dir / allele_key

        # Extract reads from ALL contigs in the cluster, then remap to
        # the peak contig.  Reads spread across multiple ladder contigs
        # due to length variation; Clair3 needs them all aligned to a
        # single reference contig to call variants effectively.
        allele_bam = _extract_and_remap_reads(
            bam_path,
            cluster_contigs,
            contig_name,
            reference_path,
            allele_dir,
            threads,
            preset=preset,
        )

        # Run Clair3 against the single-contig reference (created during
        # remapping).  Using the single contig rather than the full ladder
        # avoids Clair3 scanning 150 empty contigs.
        contig_ref = allele_dir / f"{contig_name}.fa"
        clair3_dir = allele_dir / "clair3"
        # Split thread budget across parallel alleles to avoid CPU oversubscription
        per_allele_threads = max(1, threads // len(allele_keys))
        vcf = run_clair3(
            allele_bam,
            contig_ref,
            clair3_dir,
            model_path=clair3_model,
            platform=platform,
            threads=per_allele_threads,
            settings=settings,
        )

        filtered = filter_vcf(
            vcf, contig_ref, allele_dir, min_qual=min_qual, min_dp=min_dp, **filter_options
        )
        variants = parse_vcf_genotypes(filtered)
        evidence = phase_evidence(variants)
        sample = variants[0].get("sample") if variants else None
        haplotype, genotype_status, heterozygous = length_partition_selection(variants, evidence)
        annotate_consensus_candidate(allele_info, evidence, haplotype, sample, str(filtered))
        allele_info["allele_genotype_status"] = genotype_status
        allele_info["heterozygous_sites"] = heterozygous
        allele_info["independent_haplotype_evidence"] = len(allele_keys) > 1 and haplotype == 1
        allele_info["variant_filter"] = filter_provenance
        return allele_key, filtered

    # Process both alleles in parallel when they are independent
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(_process_allele, k): k for k in allele_keys}
        results: dict[str, Path] = {}
        for future in futures:
            key, vcf_path = future.result()
            results[key] = vcf_path

    return results
