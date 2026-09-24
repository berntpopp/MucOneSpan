"""Self-contained HTML report generation.

Requires the ``jinja2`` package, which is an optional dependency.
Install with: ``pip install muc_one_span[report]``
"""

from __future__ import annotations

import math
from datetime import datetime, timezone
from pathlib import Path

try:
    from jinja2 import Environment, PackageLoader

    _HAS_JINJA2 = True
except ImportError:
    _HAS_JINJA2 = False

from typing import Any

from muc_one_span.clinical_gates import mutation_supported as _mutation_supported
from muc_one_span.nomenclature import enrich_mutation_record
from muc_one_span.report_assets import (
    REPORT_IGV_MODES,
    REPORT_IGV_OFF,
    igv_payload,
    igv_provenance,
)
from muc_one_span.version import __version__


def _enrich_mutation_nomenclature(mutation: dict) -> dict:
    """Enrich a detected mutation dict with HGVS cDNA and repeat form."""
    return enrich_mutation_record(mutation)


def _exact_match_percentage(value: Any) -> float | None:
    """Validate the producer's 0-100 percent contract without guessing units."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value) if 0 <= value <= 100 and math.isfinite(value) else None


def _execution_context(summary: dict, execution_status: dict | None) -> dict[str, Any]:
    """Prefer explicitly supplied provenance; never guess a sidecar location."""
    record = execution_status if execution_status is not None else summary.get("run_status")
    if record is None:
        return {"label": "Execution status unavailable (legacy input)", "warning": None}
    status = record.get("status") if isinstance(record, dict) else None
    if not isinstance(status, str):
        status = None
    labels = {
        "completed": "Completed",
        "analysis_completed": "Analysis completed; report generated",
        "execution_failed": "Execution failed",
        "interrupted": "Interrupted",
        "insufficient_evidence": "Insufficient evidence",
        "running": "Running",
    }
    label = (
        labels.get(status, "Unknown execution status")
        if status is not None
        else "Unknown execution status"
    )
    warning = None
    if status not in {"completed", "analysis_completed"}:
        warning = (
            f"{label}. Execution does not establish a completed analysis. "
            "Existing mutation evidence is retained; absence of detected mutations "
            "cannot establish a negative result."
        )
    return {"label": label, "warning": warning}


def compute_clinical_decision(
    summary: dict[str, Any], *, execution_status: dict | None = None
) -> dict[str, Any]:
    """Derive 3-state clinical decision support banner and multiplicity caveats."""
    execution = _execution_context(summary, execution_status)
    classifications = summary.get("classifications", {})
    alleles = summary.get("alleles", {})

    pathogenic_mutations: list[dict[str, Any]] = []
    uncertain_mutations: list[dict[str, Any]] = []
    for allele_key, acls in classifications.items():
        if isinstance(acls, dict):
            for mut in acls.get("mutations", []):
                if isinstance(mut, dict):
                    mut_copy = dict(mut)
                    mut_copy["allele"] = allele_key
                    frameshift = mut.get("frameshift") is True
                    loc_ok = mut.get("localization_status") != "ambiguous"
                    supp_ok = _mutation_supported(mut)
                    if frameshift and loc_ok and supp_ok:
                        pathogenic_mutations.append(mut_copy)
                    else:
                        uncertain_mutations.append(mut_copy)

    a1 = alleles.get("allele_1", {}) if isinstance(alleles, dict) else {}
    a2 = alleles.get("allele_2", {}) if isinstance(alleles, dict) else {}
    total_reads = (a1.get("reads", 0) or 0) + (a2.get("reads", 0) or 0)
    low_coverage = total_reads < 30 and (bool(a1) or bool(a2))

    ambiguous_bases = sum(
        acls.get("ambiguous_bases", 0)
        for acls in classifications.values()
        if isinstance(acls, dict)
    )

    reconstruction_reasons: list[str] = []
    if not bool(alleles.get("homozygous")):
        for a_key, a_info in (("Allele 1", a1), ("Allele 2", a2)):
            if not a_info:
                continue
            if a_info.get("independent_haplotype_evidence") is False:
                reconstruction_reasons.append(
                    f"{a_key}: Reconstruction incomplete; independent biological haplotype evidence not established."
                )
            elif a_info.get("phase_status") in (
                "unphased",
                "missing_phase_set",
                "disconnected_phase_sets",
                "conflicting_variant_records",
            ):
                reconstruction_reasons.append(
                    f"{a_key}: Genotype phase is unphased or conflicting."
                )
            elif a_info.get("reconstruction_status") in (
                "not_separately_resolved",
                "candidate_reference_confidence_unverified",
            ) and a_info.get("candidate_duplicate_of"):
                reconstruction_reasons.append(
                    f"{a_key}: Candidate reconstruction not separately resolved."
                )

    if pathogenic_mutations:
        state = "PATHOGENIC"
        title = "Pathogenic Variant Detected (ADTKD-MUC1)"
        badge_label = "PATHOGENIC"
        badge_class = "badge-danger"
        banner_class = "decision-pathogenic"
        details: list[str] = []
        for m in pathogenic_mutations:
            m_name = m.get("name") or m.get("mutation_name", "Unknown variant")
            rep_idx = m.get("repeat_index", "N/A")
            h_form = m.get("hgvs_repeat_form") or m.get("repeat_relative_coordinate", "")
            c_form = m.get("hgvs_cdna", "")
            supp = m.get("support_status", "")
            allele_name = str(m.get("allele", "allele")).replace("_", " ").title()
            det = f"{allele_name}: {m_name} at repeat unit {rep_idx}"
            if h_form:
                det += f" ({h_form})"
            if c_form and c_form != "transcript_coordinate_unresolved":
                det += f" [{c_form}]"
            if supp:
                det += f" - Evidence: {supp}"
            details.append(det)
        for m in uncertain_mutations:
            m_name = m.get("name") or m.get("mutation_name", "Unknown variant")
            rep_idx = m.get("repeat_index", "N/A")
            allele_name = str(m.get("allele", "allele")).replace("_", " ").title()
            details.append(
                f"{allele_name}: Additional uncertain variant ({m_name} at repeat {rep_idx}) observed."
            )
        summary_text = (
            "A pathogenic frameshift variant was identified in the MUC1 VNTR region. "
            "This finding is consistent with autosomal dominant tubulointerstitial "
            "kidney disease (ADTKD-MUC1)."
        )
        recommendations = (
            "Recommend genetic counseling and nephrology clinical correlation. "
            "Cascade variant testing is available for at-risk family members."
        )
    elif (
        low_coverage
        or ambiguous_bases > 10
        or execution["warning"] is not None
        or bool(uncertain_mutations)
        or bool(reconstruction_reasons)
    ):
        state = "INCONCLUSIVE"
        title = "Inconclusive / Quality Warning"
        badge_label = "INCONCLUSIVE"
        badge_class = "badge-warning"
        banner_class = "decision-inconclusive"
        reasons: list[str] = []
        if execution["warning"]:
            reasons.append(execution["warning"])
        if low_coverage:
            reasons.append(
                f"Total read depth ({total_reads} reads) is below diagnostic threshold (30 reads)."
            )
        if ambiguous_bases > 10:
            reasons.append(
                f"High number of ambiguous consensus bases ({ambiguous_bases}) detected."
            )
        for m in uncertain_mutations:
            m_name = m.get("name") or m.get("mutation_name", "Unknown variant")
            rep_idx = m.get("repeat_index", "N/A")
            allele_name = str(m.get("allele", "allele")).replace("_", " ").title()
            reasons.append(
                f"{allele_name}: Observed sequence variant ({m_name} at repeat {rep_idx}) "
                "is inconclusive (in-frame or ambiguous localization/support)."
            )
        if reconstruction_reasons:
            reasons.extend(reconstruction_reasons)
        if not reasons:
            reasons.append("Quality control metrics did not meet validation standards.")
        summary_text = (
            "The test result is inconclusive due to execution, quality or coverage limitations. "
            "No definitive clinical call can be rendered."
        )
        details = reasons
        recommendations = (
            "Orthogonal diagnostic validation or repeat sequencing with higher target "
            "depth is recommended prior to clinical decision-making."
        )
    else:
        state = "NO_PATHOGENIC_VARIANT_DETECTED"
        title = "No Pathogenic Variant Detected"
        badge_label = "NEGATIVE"
        badge_class = "badge-success"
        banner_class = "decision-negative"
        summary_text = (
            "No known ADTKD-MUC1 pathogenic variants (dupC, dupA, or related frameshifts) "
            "were detected across the reconstructed MUC1 VNTR alleles."
        )
        details = [
            f"Allele 1: {a1.get('length', 'N/A')} repeats ({a1.get('canonical_repeats', 'N/A')} canonical units) - {a1.get('reads', 0)} reads",
            f"Allele 2: {a2.get('length', 'N/A')} repeats ({a2.get('canonical_repeats', 'N/A')} canonical units) - {a2.get('reads', 0)} reads",
        ]
        recommendations = (
            "A negative result significantly reduces the likelihood of ADTKD-MUC1 caused by "
            "VNTR frameshift mutations. It does not exclude variants outside the VNTR or other "
            "genetic causes of kidney disease."
        )

    multiplicity_caveat = None
    l1 = a1.get("length")
    l2 = a2.get("length")
    is_single_length = (
        (l1 is not None and l2 is not None and l1 == l2)
        or bool(alleles.get("homozygous"))
        or (bool(a1) and not bool(a2))
    )
    if is_single_length and not pathogenic_mutations:
        multiplicity_caveat = (
            "Single allele length observed; second allele not established. "
            "Preferential PCR amplification or allelic drop-out cannot be excluded; "
            "apparent homozygosity should be interpreted with clinical caution."
        )
    elif is_single_length and pathogenic_mutations:
        multiplicity_caveat = (
            "Single allele length observed with pathogenic mutation; second allele not established. "
            "A second unamplified or co-migrating wild-type or mutated allele cannot be excluded."
        )

    return {
        "execution": execution,
        "state": state,
        "title": title,
        "badge_label": badge_label,
        "badge_class": badge_class,
        "banner_class": banner_class,
        "summary": summary_text,
        "details": details,
        "recommendations": recommendations,
        "multiplicity_caveat": multiplicity_caveat,
    }


def generate_report(
    summary: dict,
    output_path: Path,
    sample_name: str = "unknown",
    tool_versions: dict[str, str] | None = None,
    detailed_repeats: dict | None = None,
    *,
    report_igv: str = REPORT_IGV_OFF,
    bam_path: Path | None = None,
    vcf_path: Path | None = None,
    fasta_path: Path | None = None,
    bed_path: Path | None = None,
    vcf_paths: dict[str, Path] | None = None,
    execution_status: dict | None = None,
) -> Path:
    """Render pipeline results as a self-contained HTML report.

    Args:
        summary: Pipeline result dictionary containing ``alleles``,
            ``classifications``, ``tool_versions``, and optionally
            ``pipeline_version``.
        output_path: Destination path for the HTML file.  Parent
            directories are created automatically.
        sample_name: Human-readable sample identifier shown in the
            report header.
        tool_versions: Optional mapping of tool name to version string.
            Falls back to ``summary["tool_versions"]`` when not given.
        detailed_repeats: Optional per-allele repeat classification
            details (``repeats``, ``mutations_detected``, ``confidence``).
            When provided, a collapsible "Detailed Repeat Table" section
            is included.
        report_igv: One of ``"embedded"``, ``"sidecar"``, or ``"off"``.
        bam_path: Optional BAM track for alignment visualization.
        vcf_path: Optional VCF track for variant visualization.
        fasta_path: Optional reference FASTA for alignment visualization.
        bed_path: Optional BED region for alignment visualization.
        vcf_paths: Allele-keyed VCF tracks, sorted by key and deduplicated by path.
            When supplied (including empty), takes precedence over vcf_path.
        execution_status: Authoritative execution provenance, overriding summary
            run_status. Missing legacy status retains prior evidence decisions.
            analysis_completed is reserved for the successful pipeline analysis
            before its report callback has returned; it is not terminal status.

    Returns:
        The resolved *output_path* after writing the report.

    Raises:
        ImportError: If Jinja2 is not installed.
    """
    if not _HAS_JINJA2:
        raise ImportError(
            "Jinja2 is required for report generation. "
            "Install with: pip install muc_one_span[report]"
        )

    if report_igv not in REPORT_IGV_MODES:
        raise ValueError(
            f"Invalid report_igv mode {report_igv!r}; expected one of {REPORT_IGV_MODES}"
        )

    env = Environment(
        loader=PackageLoader("muc_one_span", "templates"),
        autoescape=True,
    )
    env.filters["exact_match_percentage"] = _exact_match_percentage
    template = env.get_template("report.html.j2")

    versions = tool_versions or summary.get("tool_versions", {})

    # Enrich classifications with HGVS nomenclature
    classifications = summary.get("classifications", {})
    enriched_classifications = {}
    for akey, acls in classifications.items():
        if isinstance(acls, dict):
            cls_copy = dict(acls)
            if "mutations" in cls_copy:
                cls_copy["mutations"] = [
                    _enrich_mutation_nomenclature(m) for m in cls_copy["mutations"]
                ]
            enriched_classifications[akey] = cls_copy
        else:
            enriched_classifications[akey] = acls

    summary_copy = dict(summary)
    summary_copy["classifications"] = enriched_classifications

    # Build IGV context if requested and reference FASTA is available
    igv_context = None
    igv_payload_b64 = None
    igv_provenance_str = None
    if report_igv != REPORT_IGV_OFF and fasta_path is not None:
        from muc_one_span.report_igv import build_igv_context

        contig_names = [
            allele["contig_name"]
            for _, allele in sorted(summary.get("alleles", {}).items())
            if isinstance(allele, dict) and allele.get("contig_name")
        ]
        igv_context = build_igv_context(
            output_dir=output_path.parent,
            fasta_path=fasta_path,
            bam_path=bam_path,
            vcf_path=vcf_path,
            vcf_paths=vcf_paths,
            bed_path=bed_path,
            report_igv=report_igv,
            flanking=0,
            contig_names=contig_names or None,
        )
        igv_payload_b64 = igv_payload(report_igv)
        igv_provenance_str = igv_provenance(report_igv)

    decision = compute_clinical_decision(summary_copy, execution_status=execution_status)

    html = template.render(
        sample_name=sample_name,
        summary=summary_copy,
        clinical_decision=decision,
        multiplicity_caveat=decision.get("multiplicity_caveat"),
        detailed_repeats=detailed_repeats,
        tool_versions=versions,
        pipeline_version=summary.get("pipeline_version", __version__),
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
        igv=igv_context,
        igv_payload=igv_payload_b64,
        igv_provenance=igv_provenance_str,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    return output_path
