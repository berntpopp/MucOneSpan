"""IGV Desktop session: relative paths, every track and every locus row as a region."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from muc_one_span.hybrid.igv_session import locus_string, write_igv_session


def test_session_lists_tracks_in_order_and_loci_as_regions(tmp_path: Path) -> None:
    names = ["igv_reference.fa", "genes.bed", "mutations.bed", "mapping.bam", "loci.bed"]
    fasta, genes, mutations, bam, loci = (tmp_path / n for n in names)
    loci.write_text("c\t10\t20\tMUC1_VNTR_allele_1\nc\t14\t15\tallele_1:repeat_3:X:dupC\n")
    path = write_igv_session(
        tmp_path,
        fasta,
        [("MUC1 gene", genes), ("Detected mutations", mutations)],
        bam,
        loci,
        locus_string("c", 14, 15, 60),
    )
    root = ET.parse(path).getroot()
    assert root.get("genome") == "igv_reference.fa" and root.get("locus") == "c:1-75"
    resources = [(r.get("name"), r.get("path")) for r in root.iter("Resource")]
    assert resources == [
        ("MUC1 gene", "genes.bed"),
        ("Detected mutations", "mutations.bed"),
        ("Assigned reads", "mapping.bam"),
    ]
    regions = [
        (r.get("start"), r.get("end"), r.get("description")) for r in root.iter("RegionOfInterest")
    ]
    assert regions == [("10", "20", "MUC1_VNTR_allele_1"), ("14", "15", "allele_1:repeat_3:X:dupC")]


def test_locus_string_is_one_based_with_context() -> None:
    assert locus_string("c", 99, 100, 10) == "c:90-110"
