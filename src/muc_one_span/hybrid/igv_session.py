"""IGV Desktop session for the hybrid IGV report's display reference and tracks.

``igv/igv_session.xml`` opens the same view as the HTML report in IGV Desktop: the
display reference (``igv_reference.fa`` as the genome), the annotation tracks, the
assigned-read alignment and every allele and mutation as a region of interest
(Regions > Region Navigator). Paths are relative to the session file, so the
``igv/`` folder can be moved or shared as a whole.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

SESSION_FILE = "igv_session.xml"
SESSION_VERSION = "8"
ALIGNMENT_TRACK = "Assigned reads"


def locus_string(contig: str, start: int, end: int, context_bp: int) -> str:
    """1-based IGV locus for a 0-based half-open interval plus ``context_bp`` each side."""
    return f"{contig}:{max(start + 1 - context_bp, 1)}-{end + context_bp}"


def write_igv_session(
    directory: Path,
    fasta: Path,
    tracks: list[tuple[str, Path]],
    bam: Path,
    loci_bed: Path,
    initial_locus: str,
) -> Path:
    """Write ``igv_session.xml`` (relative paths) into ``directory``; return its path."""
    session = ET.Element(
        "Session",
        genome=fasta.relative_to(directory).as_posix(),
        locus=initial_locus,
        version=SESSION_VERSION,
    )
    resources = ET.SubElement(session, "Resources")
    for name, path in [*tracks, (ALIGNMENT_TRACK, bam)]:
        ET.SubElement(resources, "Resource", name=name, path=path.relative_to(directory).as_posix())
    regions = ET.SubElement(session, "Regions")
    for line in loci_bed.read_text(encoding="utf-8").splitlines():
        contig, start, end, name = line.split("\t")[:4]
        ET.SubElement(
            regions,
            "RegionOfInterest",
            chromosome=contig,
            start=start,
            end=end,
            description=name,
        )
    ET.indent(session)
    path = directory / SESSION_FILE
    ET.ElementTree(session).write(path, encoding="UTF-8", xml_declaration=True)
    return path
