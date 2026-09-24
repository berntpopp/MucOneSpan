"""Interactive IGV alignment report generation (embedded or sidecar)."""

from __future__ import annotations

import base64
import gzip
import json
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


def _reference_lengths(fasta_path: Path) -> dict[str, int]:
    """Read reference contig names and lengths from its index or sequence."""
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

    if not contig_lengths:
        raise ValueError(f"Reference FASTA has no contigs: {fasta_path}")
    return contig_lengths


def _selected_contigs(contig_lengths: dict[str, int], contig_names: list[str] | None) -> list[str]:
    """Validate requested contigs, retaining stable order and removing duplicates."""
    selected_contigs = list(dict.fromkeys(contig_names or list(contig_lengths)[:1]))
    missing = [c for c in selected_contigs if c not in contig_lengths]
    if missing:
        raise ValueError(f"Requested report contigs absent from reference: {', '.join(missing)}")

    return selected_contigs


def create_locus_bed(
    fasta_path: Path,
    bed_path: Path,
    contig_names: list[str] | None = None,
    locus_name: str = "MUC1_VNTR",
) -> Path:
    """Create a 1-interval BED file covering the VNTR region of reference contigs."""
    contig_lengths = _reference_lengths(fasta_path)
    selected_contigs = _selected_contigs(contig_lengths, contig_names)

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


def _validate_vcf_payloads(output_path: Path, configs: list[dict[str, str]]) -> None:
    """Reject create_report versions that concatenate VCF records into headers."""
    headers = {}
    for config in configs:
        if config["type"] != "variant":
            continue
        path = Path(config["url"])
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#CHROM\t"):
                    headers[config["name"]] = line.rstrip("\r\n")
                    break
    if not headers:
        return
    _, _, dictionary = extract_igv_fragments(output_path.read_text(encoding="utf-8"))
    sessions = json.loads(dictionary.strip().rstrip(";"))
    for session_uri in sessions.values():
        session = json.loads(_decode_data_uri(session_uri))
        for track in session.get("tracks", []):
            name = track.get("name")
            if track.get("type") != "variant" or name not in headers:
                continue
            payload = _decode_data_uri(track["url"])
            embedded_header = next(
                (line for line in payload.splitlines() if line.startswith("#CHROM\t")), None
            )
            if embedded_header != headers[name]:
                output_path.unlink(missing_ok=True)
                raise ValueError(
                    f"Malformed VCF track {name!r} from create_report: column header changed. "
                    "Use a working igv-reports version (1.13.0 verified); 1.16.0 can "
                    "concatenate the first variant into its header."
                )


def _decode_data_uri(uri: str) -> str:
    """Decode igv-reports inline JSON/VCF data, supporting gzip and plain base64."""
    metadata, encoded = uri.split(",", 1)
    data = base64.b64decode(encoded)
    if "gzip" in metadata:
        data = gzip.decompress(data)
    return data.decode("utf-8")


def run_igv_report(
    bed_file: str | Path,
    fasta_file: str | Path,
    output_html: str | Path,
    bam_file: str | Path | None = None,
    vcf_file: str | Path | None = None,
    flanking: int = 0,
    *,
    report_igv: str = DEFAULT_REPORT_IGV,
    vcf_paths: dict[str, Path] | None = None,
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
        vcf_paths: Named VCFs; when supplied overrides singular vcf_file.
            Shared paths produce one explicitly shared track.

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
    requested_vcfs = (
        vcf_paths
        if vcf_paths is not None
        else ({"variants": Path(vcf_file)} if vcf_file is not None else {})
    )
    grouped: dict[Path, list[str]] = {}
    for allele, path in sorted(requested_vcfs.items()):
        resolved = Path(path).resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"Requested VCF track does not exist: {path}")
        grouped.setdefault(resolved, []).append(allele.replace("_", " ").title())
    configs = []
    for path, labels in grouped.items():
        name = (
            ("Variants" if labels == ["Variants"] else f"{labels[0]} variants")
            if len(labels) == 1
            else f"Shared variants ({', '.join(labels)})"
        )
        configs.append({"url": str(path), "name": name, "format": "vcf", "type": "variant"})
    if bam_file is not None:
        if not Path(bam_file).is_file():
            raise FileNotFoundError(f"Requested BAM track does not exist: {bam_file}")
        configs.append(
            {
                "url": str(Path(bam_file).resolve()),
                "name": "Alignments",
                "format": "bam",
                "type": "alignment",
            }
        )
    cmd.extend(["--output", str(output_path)])

    logger.info("Executing create_report for %s", output_html)
    with tempfile.TemporaryDirectory(prefix="muconespan-igv-tracks-") as track_dir:
        if configs:
            config_path = Path(track_dir) / "tracks.json"
            config_path.write_text(json.dumps(configs), encoding="utf-8")
            cmd.extend(["--track-config", str(config_path)])
        run_tool(cmd)
    _validate_vcf_payloads(output_path, configs)

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


def preflight_igv_report(report_igv: str, work_dir: Path) -> None:
    """Fail before analysis when the requested IGV report cannot be produced.

    Runs ``create_report`` once on a synthetic one-record locus in a temporary
    directory, reusing :func:`run_igv_report` validation. A missing executable or
    a version that corrupts VCF tracks (observed with igv-reports 1.16.x) then
    stops the run before mapping instead of after the analysis has completed.
    """
    if report_igv == REPORT_IGV_OFF:
        return
    work_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="igv-preflight-", dir=work_dir) as tmp:
        root = Path(tmp)
        sequence = "ACGT" * 300
        (root / "probe.fa").write_text(f">probe\n{sequence}\n", encoding="utf-8")
        (root / "probe.bed").write_text("probe\t100\t200\tprobe\n", encoding="utf-8")
        (root / "probe.vcf").write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=probe,length=1200>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"probe\t150\t.\t{sequence[149]}\tT\t30\tPASS\t.\n",
            encoding="utf-8",
        )
        try:
            run_igv_report(
                root / "probe.bed",
                root / "probe.fa",
                root / "probe.html",
                report_igv=report_igv,
                vcf_paths={"probe": root / "probe.vcf"},
            )
        except (OSError, RuntimeError, ValueError) as exc:
            raise RuntimeError(
                f"IGV report preflight failed for --report-igv {report_igv}: {exc}. "
                "Use a working igv-reports create_report (1.13.0 verified) or --report-igv off."
            ) from exc


def build_igv_context(
    output_dir: Path,
    fasta_path: Path,
    bam_path: Path | None = None,
    vcf_path: Path | None = None,
    bed_path: Path | None = None,
    report_igv: str = DEFAULT_REPORT_IGV,
    flanking: int = 0,
    contig_names: list[str] | None = None,
    *,
    vcf_paths: dict[str, Path] | None = None,
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

    if not fasta_path.is_file():
        raise FileNotFoundError(f"Requested reference FASTA does not exist: {fasta_path}")
    if bed_path is not None and not bed_path.is_file():
        raise FileNotFoundError(f"Requested BED does not exist: {bed_path}")
    if bed_path is not None and contig_names:
        _selected_contigs(_reference_lengths(fasta_path), contig_names)
    temp_dir = None
    try:
        if bed_path is None:
            temp_dir = tempfile.TemporaryDirectory(prefix="muconespan-igv-")
            bed_path = create_locus_bed(
                fasta_path,
                Path(temp_dir.name) / "locus.bed",
                contig_names=contig_names,
            )
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
                vcf_paths=vcf_paths,
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
                vcf_paths=vcf_paths,
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
