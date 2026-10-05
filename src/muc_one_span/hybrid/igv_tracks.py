"""IGV display reference, locus rows and detected-mutation track for hybrid reports.

The hybrid engine writes no alignment of its own. When an IGV report is requested,
each allele is shown against a *display reference*: its polished consensus with every
mutated repeat unit restored to the unit's canonical parent sequence (the closest
dictionary repeat, or the template's parent unit), wrapped in the same ladder flanks
that read assignment uses (:func:`muc_one_span.hybrid.assign.hybrid_references`).
Reads carrying a detected mutation therefore show it against the reference as an
insertion, deletion or mismatch, exactly where the mutation track marks it.
Unmutated units keep the consensus sequence, so only the called events stand out.

Coordinates are 0-based half-open (BED). The changed bases come from the
unit-to-observed differences (:func:`muc_one_span.repeat_alignment.characterize_differences`),
whose 1-based positions are in the canonical unit, i.e. directly in the display
reference: a substitution or deletion marks the reference base(s) concerned, an
insertion the two reference bases around its insertion point. Insertions and
deletions are left-normalized within a homopolymer run (as aligners place them), so
the mark sits where the reads show the event. Each mutation's navigation row is its
1-bp *site*; igv-reports then sorts the reads by base there (carriers first) and
adds the requested context on each side.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.igv_gene import GeneModel, gene_bed_line, parse_region
from muc_one_span.repeat_alignment import characterize_differences

REFERENCE_PREFIX = "hybrid_"
ALLELE_LOCUS_SUFFIX = "MUC1_VNTR"


@dataclass(frozen=True)
class MutationFeature:
    """One detected mutation on a display-reference contig (0-based, half-open)."""

    contig: str
    unit_start: int
    unit_end: int
    start: int
    end: int
    site: int
    name: str

    def bed_line(self) -> str:
        """Track line: the changed reference bases."""
        return f"{self.contig}\t{self.start}\t{self.end}\t{self.name}\n"

    def locus_line(self) -> str:
        """Navigation row: the 1-bp site (reads are sorted by base there)."""
        return f"{self.contig}\t{self.site}\t{self.site + 1}\t{self.name}\n"


def _left_normalized(index: int, base: str, unit: str) -> int:
    """Move ``index`` left while the preceding unit base equals ``base``."""
    while index > 0 and unit[index - 1] == base:
        index -= 1
    return index


def variant_site(differences: list[dict[str, Any]], unit: str) -> tuple[int, int, int]:
    """Changed unit interval (0-based, half-open) and the sort site of a mutation.

    Positions are 1-based in the canonical ``unit``. A substitution marks its base. A
    deletion marks the deleted base, and an insertion the two bases around its
    insertion point, both left-normalized within a run of the same base. The site is
    the leftmost marked event: a substitution's or deleted base, or the base right
    after an insertion point (where an aligner records the insertion). An empty
    list marks the whole unit, with the site at its start.
    """
    lo, hi, site = len(unit), 0, None
    for diff in sorted(differences, key=lambda d: d["pos"]):
        p = diff["pos"]
        if diff["type"] == "insertion":
            point = _left_normalized(p - 1, diff["alt"][:1], unit)
            first, last, here = point - 1, point + 1, point
        elif diff["type"] == "deletion":
            base = _left_normalized(p - 1, diff["ref"][:1], unit)
            first, last, here = base, base + 1, base
        else:
            first, last, here = p - 1, p, p - 1
        lo, hi = min(lo, first), max(hi, last)
        site = here if site is None else min(site, here)
    lo, hi = max(lo, 0), min(hi, len(unit))
    if site is None or lo >= hi:
        return 0, len(unit), 0
    return lo, hi, min(max(site, 0), len(unit) - 1)


def _label(mutation: dict[str, Any], differences: list[dict[str, Any]]) -> str:
    if mutation.get("mutation_name"):
        return str(mutation["mutation_name"])
    parts = []
    for d in differences:
        if d["type"] == "insertion":
            parts.append(f"ins{d['pos']}{d['alt']}")
        elif d["type"] == "deletion":
            parts.append(f"del{d['pos']}{d['ref']}")
        else:
            parts.append(f"sub{d['pos']}{d['ref']}/{d['alt']}")
    return ",".join(parts)


# Repeat-unit track colours: the light-theme --clr-* borders of report.css, so the IGV
# track matches the report's repeat-structure legend.
UNIT_COLORS = {
    "pre": "29,78,216",
    "canonical": "107,114,128",
    "variant": "4,120,87",
    "after": "180,83,9",
    "mutation": "185,28,28",
}
VARIANT_PREFIX = "?"
MUTATED_SUFFIX = "m"


@dataclass(frozen=True)
class UnitFeature:
    """One classified repeat unit on a display contig (0-based, half-open)."""

    start: int
    end: int
    label: str
    category: str

    def bed_line(self, contig: str, offset: int, index: int) -> str:
        """BED9 line coloured by category (``itemRgb``)."""
        start, end = self.start + offset, self.end + offset
        color = UNIT_COLORS[self.category]
        return f"{contig}\t{start}\t{end}\t{index}:{self.label}\t0\t+\t{start}\t{end}\t{color}\n"


def unit_category(label: str, rd: RepeatDictionary) -> str:
    """Report legend category of one structure token (as the report's repeat map)."""
    if ":" in label or label.endswith(MUTATED_SUFFIX):
        return "mutation"
    if label.startswith(VARIANT_PREFIX):
        return "variant"
    if label in rd.pre_repeat_ids:
        return "pre"
    if label in rd.after_repeat_ids:
        return "after"
    return "canonical"


@dataclass(frozen=True)
class DisplayAllele:
    """Display sequence of one allele (no flank) with its mutation and unit features."""

    sequence: str
    mutations: list[MutationFeature]
    units: list[UnitFeature]


def display_allele(
    allele: str, consensus: str, classification: dict[str, Any], rd: RepeatDictionary
) -> DisplayAllele:
    """Consensus with mutated units restored, its mutation and unit features (no flank)."""
    contig = f"{REFERENCE_PREFIX}{allele}"
    by_index: dict[int, list[dict[str, Any]]] = {}
    for mutation in classification.get("mutations_detected", []):
        by_index.setdefault(mutation["repeat_index"], []).append(mutation)
    units = sorted(classification.get("repeats", []), key=lambda r: r["start"])
    missing = set(by_index) - {r["index"] for r in units}
    if missing:
        raise ValueError(f"{allele}: mutation at repeat {min(missing)} has no classified window")
    labels = classification.get("structure", "").split()
    parts: list[str] = []
    features: list[MutationFeature] = []
    unit_features: list[UnitFeature] = []
    length = previous_end = 0
    for position, unit in enumerate(units):
        label = labels[position] if position < len(labels) else str(unit.get("type", ""))
        unit_start = length + len(consensus[previous_end : unit["start"]])
        gap = consensus[previous_end : unit["start"]]
        parts.append(gap)
        length += len(gap)
        observed = consensus[unit["start"] : unit["end"]]
        previous_end = unit["end"]
        mutations = by_index.get(unit["index"])
        if not mutations:
            parts.append(observed)
            length += len(observed)
            unit_features.append(UnitFeature(unit_start, length, label, unit_category(label, rd)))
            continue
        parent = rd.repeats[mutations[0]["closest_type"]]
        for mutation in mutations:
            differences = mutation.get("differences")
            if differences is None:
                differences = characterize_differences(parent, observed)
            lo, hi, site = variant_site(differences, parent)
            name = f"{allele}:repeat_{unit['index']}:{mutation['closest_type']}:"
            features.append(
                MutationFeature(
                    contig,
                    length,
                    length + len(parent),
                    length + lo,
                    length + hi,
                    length + site,
                    name + _label(mutation, differences),
                )
            )
        parts.append(parent)
        length += len(parent)
        unit_features.append(UnitFeature(unit_start, length, label, "mutation"))
    parts.append(consensus[previous_end:])
    return DisplayAllele("".join(parts), features, unit_features)


def _shifted(feature: MutationFeature, offset: int) -> MutationFeature:
    return MutationFeature(
        feature.contig,
        feature.unit_start + offset,
        feature.unit_end + offset,
        feature.start + offset,
        feature.end + offset,
        feature.site + offset,
        feature.name,
    )


@dataclass(frozen=True)
class IgvFiles:
    """Display reference, navigation rows and annotation tracks (name, BED path)."""

    fasta: Path
    loci: Path
    tracks: list[tuple[str, Path]]
    mutation_names: list[str]
    contig_lengths: dict[str, int]


GENE_TRACK = "MUC1 gene"
UNIT_TRACK = "Repeat units"
MUTATION_TRACK = "Detected mutations"


def write_igv_inputs(
    directory: Path,
    consensus: Mapping[str, str],
    classifications: Mapping[str, dict[str, Any]],
    rd: RepeatDictionary,
    flanks: tuple[int, int],
    gene: GeneModel | None = None,
) -> IgvFiles:
    """Write the display reference, ``loci.bed`` and the annotation track BEDs.

    Each allele becomes contig ``hybrid_<allele>`` (the summary ``contig_name``):
    ``flanks[0]`` bp of left flank, the display allele and ``flanks[1]`` bp of right
    flank. Tracks: the gene model (when given), the coloured repeat units and the
    detected mutations. The locus rows (one per allele VNTR and per mutation site)
    populate the IGV navigation table.
    """
    directory.mkdir(parents=True, exist_ok=True)
    left_bp, right_bp = flanks
    left, right = rd.flanking_left[-left_bp:] if left_bp else "", rd.flanking_right[:right_bp]
    vntr = parse_region(rd.vntr_region)
    fasta_parts: list[str] = []
    loci: list[str] = []
    mutations: list[str] = []
    units: list[str] = []
    genes: list[str] = []
    lengths: dict[str, int] = {}
    names: list[str] = []
    for allele in sorted(consensus):
        shown = display_allele(allele, consensus[allele], classifications.get(allele, {}), rd)
        contig = f"{REFERENCE_PREFIX}{allele}"
        sequence = left + shown.sequence + right
        fasta_parts.append(f">{contig}\n{sequence}\n")
        lengths[contig] = len(sequence)
        end = left_bp + len(shown.sequence)
        loci.append(f"{contig}\t{left_bp}\t{end}\t{ALLELE_LOCUS_SUFFIX}_{allele}\n")
        for feature in (_shifted(f, left_bp) for f in shown.mutations):
            loci.append(feature.locus_line())
            mutations.append(feature.bed_line())
            names.append(feature.name)
        units.extend(u.bed_line(contig, left_bp, i + 1) for i, u in enumerate(shown.units))
        if gene is not None:
            genes.append(gene_bed_line(gene, vntr, left_bp, len(shown.sequence), contig))
    fasta = directory / "igv_reference.fa"
    fasta.write_text("".join(fasta_parts), encoding="utf-8")
    # itemRgb="On" makes IGV Desktop (and igv.js) colour the units by category.
    units.insert(0, f'track name="{UNIT_TRACK}" itemRgb="On"\n')
    files = {"loci.bed": loci, "genes.bed": genes, "units.bed": units, "mutations.bed": mutations}
    for name, lines in files.items():
        (directory / name).write_text("".join(lines), encoding="utf-8")
    tracks = [(UNIT_TRACK, directory / "units.bed"), (MUTATION_TRACK, directory / "mutations.bed")]
    if gene is not None:
        tracks.insert(0, (GENE_TRACK, directory / "genes.bed"))
    return IgvFiles(fasta, directory / "loci.bed", tracks, names, lengths)
