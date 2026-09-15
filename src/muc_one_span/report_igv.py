"""Interactive IGV alignment report generation (embedded or sidecar)."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import Any

from muc_one_span.report_assets import (
    DEFAULT_REPORT_IGV,
    EMPTY_SESSION_DICTIONARY,
    EMPTY_TABLE_JSON,
    IGV_REPORT_LIBRARY_MARKER,
    IGV_REPORT_TEMPLATE_PATH,
    IGV_SESSION_DICTIONARY_MARKER,
    IGV_TABLE_JSON_MARKER,
    REPORT_IGV_EMBEDDED,
    REPORT_IGV_OFF,
    REPORT_IGV_SIDECAR,
    extract_igv_fragments,
    igv_library_source,
    js_json_literal,
)
from muc_one_span.tools import run_tool

logger = logging.getLogger(__name__)


def create_locus_bed(
    fasta_path: Path,
    bed_path: Path,
    contig_names: list[str] | None = None,
    locus_name: str = "MUC1_VNTR",
) -> Path:
    """Create a 1-interval BED file covering the VNTR region of reference contigs."""
    fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
    contig_lengths: dict[str, int] = {}
    if fai_path.exists():
        for line in fai_path.read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[1].isdigit():
                contig_lengths[parts[0]] = int(parts[1])
    else:
        current_name = None
        current_len = 0
        for line in fasta_path.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                if current_name is not None:
                    contig_lengths[current_name] = current_len
                current_name = line.lstrip(">").split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_name is not None:
            contig_lengths[current_name] = current_len

    selected_contigs = []
    if contig_names:
        for c in contig_names:
            if c in contig_lengths:
                selected_contigs.append(c)
    if not selected_contigs:
        selected_contigs = list(contig_lengths.keys())[:1] if contig_lengths else ["contig_1"]

    bed_lines = []
    for c in selected_contigs:
        c_len = contig_lengths.get(c, 1600)
        # igv-reports expands regionSeq by -500. For ladder contigs (c_len > 600), start must be >= 501.
        # Ladder flanks are 500 bp, so 501 is the exact 1-based start of the VNTR sequence.
        if c_len > 600:
            vntr_start = 501
            vntr_end = max(vntr_start + 60, c_len - 500 if c_len > 1000 else c_len)
        else:
            vntr_start = 0
            vntr_end = c_len
        name = f"{locus_name}_{c}" if len(selected_contigs) > 1 else locus_name
        bed_lines.append(f"{c}\t{vntr_start}\t{vntr_end}\t{name}\n")

    bed_path.parent.mkdir(parents=True, exist_ok=True)
    bed_path.write_text("".join(bed_lines), encoding="utf-8")
    return bed_path


def run_igv_report(
    bed_file: str | Path,
    fasta_file: str | Path,
    output_html: str | Path,
    bam_file: str | Path | None = None,
    vcf_file: str | Path | None = None,
    flanking: int = 0,
    *,
    report_igv: str = DEFAULT_REPORT_IGV,
) -> Path:
    """Execute igv-reports against the controlled offline template.

    Args:
        bed_file: BED file defining the displayed region.
        fasta_file: Reference FASTA used to embed sequence data.
        output_html: Target HTML path.
        bam_file: Optional BAM alignment track.
        vcf_file: Optional VCF variant track.
        flanking: Flanking base pairs around locus.
        report_igv: One of 'embedded' or 'sidecar'.

    Returns:
        Path to the generated HTML file.
    """
    if report_igv not in (REPORT_IGV_EMBEDDED, REPORT_IGV_SIDECAR):
        raise ValueError(f"--report-igv {report_igv!r} cannot generate an IGV report.")
    if not IGV_REPORT_TEMPLATE_PATH.is_file():
        raise ValueError(f"Controlled IGV report template missing: {IGV_REPORT_TEMPLATE_PATH}")

    output_path = Path(output_html)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "create_report",
        str(bed_file),
        "--template",
        str(IGV_REPORT_TEMPLATE_PATH),
        "--flanking",
        str(flanking),
        "--fasta",
        str(fasta_file),
    ]
    tracks: list[str] = []
    for track_path in (vcf_file, bam_file):
        if track_path and Path(track_path).exists():
            tracks.append(str(track_path))
    if tracks:
        cmd.extend(["--tracks", *tracks])
    cmd.extend(["--output", str(output_path)])

    logger.info("Executing create_report for %s", output_html)
    run_tool(cmd)

    if report_igv == REPORT_IGV_SIDECAR:
        generated = output_path.read_text(encoding="utf-8")
        for marker, fallback, name in (
            (IGV_TABLE_JSON_MARKER, EMPTY_TABLE_JSON, "tableJson"),
            (IGV_SESSION_DICTIONARY_MARKER, EMPTY_SESSION_DICTIONARY, "sessionDictionary"),
        ):
            if generated.count(marker) != 1:
                raise ValueError(
                    f"Generated IGV sidecar does not contain exactly one {name} marker."
                )
            value_start = generated.index(marker) + len(marker)
            value_end = generated.find("\n", value_start)
            if value_end == -1:
                value_end = len(generated)
            fragment = generated[value_start:value_end].strip()
            safe = js_json_literal(fragment, fallback)
            generated = generated[:value_start] + safe + generated[value_end:]

        if generated.count(IGV_REPORT_LIBRARY_MARKER) != 1:
            raise ValueError("Generated IGV sidecar does not contain verified library marker.")
        source = igv_library_source().decode("utf-8")
        output_path.write_text(
            generated.replace(IGV_REPORT_LIBRARY_MARKER, source), encoding="utf-8"
        )

    return output_path


def build_igv_context(
    output_dir: Path,
    fasta_path: Path,
    bam_path: Path | None = None,
    vcf_path: Path | None = None,
    bed_path: Path | None = None,
    report_igv: str = DEFAULT_REPORT_IGV,
    flanking: int = 0,
    contig_names: list[str] | None = None,
) -> dict[str, Any]:
    """Prepare IGV report and extract context for template injection."""
    if report_igv == REPORT_IGV_OFF:
        return {
            "mode": REPORT_IGV_OFF,
            "has_igv": False,
            "igv_content": "",
            "table_json": EMPTY_TABLE_JSON,
            "session_dictionary": EMPTY_SESSION_DICTIONARY,
            "sidecar_path": None,
        }

    temp_dir = None
    if bed_path is None or not bed_path.exists():
        temp_dir = tempfile.TemporaryDirectory(prefix="muconespan-igv-")
        bed_path = create_locus_bed(
            fasta_path,
            Path(temp_dir.name) / "locus.bed",
            contig_names=contig_names,
        )

    try:
        if report_igv == REPORT_IGV_EMBEDDED:
            if temp_dir is None:
                temp_dir = tempfile.TemporaryDirectory(prefix="muconespan-igv-")
            tmp_html = Path(temp_dir.name) / "igv_report.html"
            run_igv_report(
                bed_path,
                fasta_path,
                tmp_html,
                bam_file=bam_path,
                vcf_file=vcf_path,
                flanking=flanking,
                report_igv=REPORT_IGV_EMBEDDED,
            )
            content, t_json, s_dict = extract_igv_fragments(tmp_html.read_text(encoding="utf-8"))
            return {
                "mode": REPORT_IGV_EMBEDDED,
                "has_igv": bool(content and s_dict),
                "igv_content": content,
                "table_json": js_json_literal(t_json, EMPTY_TABLE_JSON),
                "session_dictionary": js_json_literal(s_dict, EMPTY_SESSION_DICTIONARY),
                "sidecar_path": None,
            }
        else:  # REPORT_IGV_SIDECAR
            sidecar_path = output_dir / "igv_report.html"
            run_igv_report(
                bed_path,
                fasta_path,
                sidecar_path,
                bam_file=bam_path,
                vcf_file=vcf_path,
                flanking=flanking,
                report_igv=REPORT_IGV_SIDECAR,
            )
            return {
                "mode": REPORT_IGV_SIDECAR,
                "has_igv": True,
                "igv_content": "",
                "table_json": EMPTY_TABLE_JSON,
                "session_dictionary": EMPTY_SESSION_DICTIONARY,
                "sidecar_path": sidecar_path,
            }
    finally:
        if temp_dir is not None:
            temp_dir.cleanup()
