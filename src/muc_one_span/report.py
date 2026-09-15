"""Self-contained HTML report generation.

Requires the ``jinja2`` package, which is an optional dependency.
Install with: ``pip install muc_one_span[report]``
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

try:
    from jinja2 import Environment, PackageLoader

    _HAS_JINJA2 = True
except ImportError:
    _HAS_JINJA2 = False

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
    if report_igv != REPORT_IGV_OFF and fasta_path and fasta_path.exists():
        from muc_one_span.report_igv import build_igv_context

        contig_names = [
            d["contig"]
            for d in (detailed_repeats or {}).values()
            if isinstance(d, dict) and "contig" in d
        ]
        igv_context = build_igv_context(
            output_dir=output_path.parent,
            fasta_path=fasta_path,
            bam_path=bam_path,
            vcf_path=vcf_path,
            bed_path=bed_path,
            report_igv=report_igv,
            flanking=0,
            contig_names=contig_names or None,
        )
        igv_payload_b64 = igv_payload(report_igv)
        igv_provenance_str = igv_provenance(report_igv)

    html = template.render(
        sample_name=sample_name,
        summary=summary_copy,
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
