"""Literature-compatible and HGVS-aligned nomenclature for MUC1 VNTR variants.

Coordinates and orientation:
- MUC1 is located on chromosome 1 on the genomic minus (reverse) strand.
- Transcription proceeds right-to-left in genomic coordinates: the coding sequence
  runs 5' to 3' along the canonical 60 bp repeat unit.
- Positions are 1-based and inclusive in the coding orientation.
- Insertions are represented as the zero-width interbase interval between two
  adjacent bases (end == start - 1).
- Per HGVS recommendations (https://hgvs-nomenclature.org/stable/), repeated
  sequence notation (c.seq[N]) in coding sequences is restricted to repeat units
  that are multiples of 3. Frameshifting indels (e.g. 59dupC) are described as
  duplications or insertions rather than repeat copy transitions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, ClassVar

CANONICAL_UNIT: str = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"
UNIT_LENGTH: int = 60
MUC1_TRANSCRIPT_ID: str = "NM_001204286.1"
MUC1_GENOMIC_REGION_HG38: str = "chr1:155188487-155192239"

KNOWN_VARIANTS: dict[str, str] = {
    "59dupC": "Kirby et al. 2013 (PMID:23396133); Wenzel et al. 2018 (PMID:29520014)",
    "56_59dupCCCC": "Vrbacka et al. 2025 (doi:10.1101/2024.11.14.623419)",
    "58_59insG": "Olinger et al. 2020 (PMID:32647000)",
    "60dupA": "Olinger et al. 2020 (PMID:32647000)",
    "55delinsAT": "Olinger et al. 2020 (PMID:32647000)",
    "54_56delinsAT": "Olinger et al. 2020 (PMID:32647000)",
    "1_5delGCCCA": "Saei et al. 2023 (PMID:37456840)",
    "30_31insCAGGCCGGCCCCGGGCTCCGGACAC": "Saei et al. 2023 (PMID:37456840)",
    "23dupC": "Vrbacka et al. 2025 (doi:10.1101/2024.11.14.623419)",
    "57_58insG": "Vrbacka et al. 2025 (doi:10.1101/2024.11.14.623419)",
    "53_54insG": "Vrbacka et al. 2025 (doi:10.1101/2024.11.14.623419)",
    "54dupG": "Vrbacka et al. 2025 (doi:10.1101/2024.11.14.623419)",
    "53_54insA": "Vrbacka et al. 2025 (doi:10.1101/2024.11.14.623419)",
}

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(sequence: str) -> str:
    """Reverse-complement a DNA nucleotide sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def trim_edit(unit: str, start: int, end: int, inserted: str) -> tuple[int, int, str]:
    """Reduce an edit to its minimal span by trimming shared flanking bases."""
    deleted = unit[start - 1 : end]
    while deleted and inserted and deleted[0] == inserted[0]:
        deleted, inserted = deleted[1:], inserted[1:]
        start += 1
    while deleted and inserted and deleted[-1] == inserted[-1]:
        deleted, inserted = deleted[:-1], inserted[:-1]
        end -= 1
    return start, end, inserted


def normalise(unit: str, start: int, end: int, inserted: str) -> tuple[int, int, str]:
    """Shift an edit as far 3' as the sequence allows (HGVS 3'-most rule).

    In the MUC1 coding orientation, 3' corresponds to increasing coordinate
    indices. A delins is anchored and never shifted.
    """
    if end > len(unit) or start < 1:
        return start, end, inserted

    start, end, inserted = trim_edit(unit, start, end, inserted)
    length = len(unit)

    if inserted and end >= start:
        return start, end, inserted

    if not inserted:
        # Deletion: roll 3' while base exiting 5' equals base entering 3'
        while end < length and unit[start - 1] == unit[end]:
            start, end = start + 1, end + 1
        return start, end, inserted

    # Insertion: roll 3' while next reference base equals first inserted base
    while end < length and inserted[0] == unit[end]:
        inserted = inserted[1:] + inserted[0]
        start, end = start + 1, end + 1
    return start, end, inserted


def roll_5prime(unit: str, start: int, end: int, inserted: str) -> tuple[int, int, str]:
    """Shift an edit as far 5' as possible to determine the ambiguity window."""
    if end > len(unit) or start < 1:
        return start, end, inserted

    if not inserted:
        while start > 1 and unit[end - 1] == unit[start - 2]:
            start, end = start - 1, end - 1
        return start, end, inserted

    while start > 1 and inserted[-1] == unit[start - 2]:
        inserted = inserted[-1] + inserted[:-1]
        start, end = start - 1, end - 1
    return start, end, inserted


def ambiguity_interval(unit: str, start: int, end: int, inserted: str) -> tuple[int, int] | None:
    """Calculate the 1-based inclusive interval of equivalent variant placement."""
    if end > len(unit) or start < 1:
        return None

    start, end, inserted = trim_edit(unit, start, end, inserted)
    if inserted and end >= start:
        return None

    low_start, _, _ = roll_5prime(unit, start, end, inserted)
    high_start, high_end, _ = normalise(unit, start, end, inserted)

    if inserted:
        low, high = low_start, high_start - 1
    else:
        low, high = low_start, high_end

    if high <= low:
        return None
    return low, high


def _tract_at(unit: str, position: int) -> tuple[int, str, int] | None:
    """Find the homopolymer run covering a 1-based position."""
    if not 1 <= position <= len(unit):
        return None
    base = unit[position - 1]
    start = position
    while start > 1 and unit[start - 2] == base:
        start -= 1
    end = position
    while end < len(unit) and unit[end] == base:
        end += 1
    count = end - start + 1
    return (start, base, count) if count >= 2 else None


def repeat_form(unit: str, start: int, end: int, inserted: str) -> str | None:
    """Express an indel as a tract copy-number transition (e.g. 53C[7]>53C[8])."""
    if end > len(unit) or start < 1:
        return None

    start, end, inserted = trim_edit(unit, start, end, inserted)
    if inserted and end >= start:
        return None

    if inserted:
        if len(set(inserted)) != 1:
            return None
        anchor = start - 1
        delta = len(inserted)
    else:
        deleted = unit[start - 1 : end]
        if not deleted or len(set(deleted)) != 1:
            return None
        anchor = start
        delta = -(end - start + 1)

    tract = _tract_at(unit, anchor)
    if tract is None:
        return None

    tract_start, base, count = tract
    expected_base = inserted[0] if inserted else unit[start - 1]
    if base != expected_base:
        return None

    new_count = count + delta
    return f"{tract_start}{base}[{count}]>{tract_start}{base}[{new_count}]"


def classify_event(start: int, end: int, inserted: str) -> str:
    """Classify the molecular event type of a normalized edit."""
    if end < start:
        return "duplication" if is_duplication(CANONICAL_UNIT, start, inserted) else "insertion"
    if not inserted:
        return "deletion"
    if start == end and len(inserted) == 1:
        return "substitution"
    return "delins"


def is_duplication(unit: str, start: int, inserted: str) -> bool:
    """Determine whether an insertion duplicates immediately preceding bases."""
    left = start - 1
    ins_len = len(inserted)
    if left < ins_len:
        return False
    return unit[left - ins_len : left] == inserted


def name_edit(unit: str, start: int, end: int, inserted: str) -> tuple[str, str]:
    """Name an edit on the repeat unit in the coding frame.

    Returns:
        tuple[str, str]: (canonical_name, event_type).
        Examples: ('59dupC', 'duplication'), ('58_59insG', 'insertion'),
                  ('54_56delinsAT', 'delins').
    """
    was_insertion = end == start - 1
    start, end, inserted = normalise(unit, start, end, inserted)

    if was_insertion or end < start:
        left = start - 1
        if is_duplication(unit, start, inserted):
            if len(inserted) == 1:
                return f"{left}dup{inserted}", "duplication"
            dup_start = left - len(inserted) + 1
            return f"{dup_start}_{left}dup{inserted}", "duplication"
        return f"{left}_{left + 1}ins{inserted}", "insertion"

    deleted = unit[start - 1 : end]
    if not inserted:
        span = str(start) if start == end else f"{start}_{end}"
        return f"{span}del{deleted}", "deletion"

    if start == end and len(inserted) == 1:
        return f"{start}{deleted}>{inserted}", "substitution"

    span = str(start) if start == end else f"{start}_{end}"
    return f"{span}delins{inserted}", "delins"


def format_hgvs_cdna(
    canonical_name: str,
    event_type: str,
    transcript: str = MUC1_TRANSCRIPT_ID,
    *,
    transcript_offset: int | None = None,
    allow_unmapped: bool = False,
) -> str:
    """Format official HGVS cDNA representation for the variant.

    To adhere to HGVS 20.05 and prevent misleading coordinates (e.g. c.59dupC in
    the signal peptide), transcript coordinate is emitted only when an authentic
    transcript offset is provided or allow_unmapped is explicitly enabled.
    """
    if not allow_unmapped and transcript_offset is None:
        return "transcript_coordinate_unresolved"
    offset_str = f"{transcript_offset}_" if transcript_offset is not None else ""
    return f"{transcript}:c.{offset_str}{canonical_name}"


def format_repeat_relative_coordinate(
    repeat_index: int | str | None,
    canonical_name: str,
) -> str:
    """Format primary invariant clinical coordinate: repeat_{idx}:c.{edit}."""
    idx = repeat_index if repeat_index is not None else "?"
    return f"repeat_{idx}:c.{canonical_name}"


@dataclass(frozen=True)
class NomenclatureRecord:
    """Structured representation of a named MUC1-VNTR variant."""

    canonical_name: str
    event_type: str
    unit_symbol: str
    repeat_position: int | None
    ambiguity_interval: tuple[int, int] | None
    repeat_form: str | None
    hgvs_cdna: str
    confidence_tier: str
    is_known_variant: bool
    literature_citation: str | None
    repeat_relative_coordinate: str = ""
    transcript_coordinate: str = "transcript_coordinate_unresolved"

    TIER_A: ClassVar[str] = "Tier_A"
    TIER_B: ClassVar[str] = "Tier_B"
    TIER_C: ClassVar[str] = "Tier_C"

    def to_dict(self) -> dict:
        """Convert to clinical reporting dictionary."""
        return {
            "canonical_name": self.canonical_name,
            "event_type": self.event_type,
            "unit_symbol": self.unit_symbol,
            "repeat_position": self.repeat_position,
            "ambiguity_interval": list(self.ambiguity_interval)
            if self.ambiguity_interval
            else None,
            "repeat_form": self.repeat_form,
            "hgvs_cdna": self.hgvs_cdna,
            "repeat_relative_coordinate": self.repeat_relative_coordinate,
            "transcript_coordinate": self.transcript_coordinate,
            "confidence_tier": self.confidence_tier,
            "is_known_variant": self.is_known_variant,
            "literature_citation": self.literature_citation,
        }


def name_variant_call(
    start: int,
    end: int,
    inserted: str,
    unit_symbol: str = "X",
    unit_sequence: str = CANONICAL_UNIT,
    support_reads: int = 10,
    independent_sources: int = 1,
    repeat_index: int | str | None = None,
) -> NomenclatureRecord:
    """Generate a validated NomenclatureRecord for a detected variant."""
    canon_name, event_type = name_edit(unit_sequence, start, end, inserted)
    interval = ambiguity_interval(unit_sequence, start, end, inserted)
    rep_form = repeat_form(unit_sequence, start, end, inserted)
    hgvs = format_hgvs_cdna(canon_name, event_type)
    citation = KNOWN_VARIANTS.get(canon_name)
    is_known = citation is not None
    repeat_rel = format_repeat_relative_coordinate(repeat_index, canon_name)

    if is_known and support_reads >= 5 and independent_sources >= 1:
        tier = NomenclatureRecord.TIER_A
    elif support_reads >= 3:
        tier = NomenclatureRecord.TIER_B
    else:
        tier = NomenclatureRecord.TIER_C

    # Extract primary numeric position for indexing
    pos = start if start <= len(unit_sequence) else None

    return NomenclatureRecord(
        canonical_name=canon_name,
        event_type=event_type,
        unit_symbol=unit_symbol,
        repeat_position=pos,
        ambiguity_interval=interval,
        repeat_form=rep_form,
        hgvs_cdna=hgvs,
        confidence_tier=tier,
        is_known_variant=is_known,
        literature_citation=citation,
        repeat_relative_coordinate=repeat_rel,
        transcript_coordinate="transcript_coordinate_unresolved",
    )


def enrich_mutation_record(
    mutation: dict[str, Any],
    repeat_dict: Any | None = None,
) -> dict[str, Any]:
    """Enrich a detected mutation dict with HGVS cDNA and repeat form."""
    enriched = dict(mutation)
    mut_name = str(enriched.get("mutation_name", ""))
    closest_type = str(enriched.get("closest_type", "X"))
    vcf_support = bool(enriched.get("vcf_support", False))
    support_reads = 10 if vcf_support else (3 if enriched.get("frameshift") else 1)
    repeat_idx = enriched.get("repeat_index")

    start: int | None = None
    end: int | None = None
    inserted: str | None = None

    # 1. Resolve from repeat dictionary definitions if available
    if repeat_dict is None:
        try:
            from muc_one_span.config import _bundled_repeats_path, load_repeat_dictionary

            repeat_dict = load_repeat_dictionary(_bundled_repeats_path())
        except Exception:
            repeat_dict = None

    if repeat_dict and hasattr(repeat_dict, "mutations") and mut_name in repeat_dict.mutations:
        mdef = repeat_dict.mutations[mut_name]
        changes = mdef.get("changes", [])
        if changes:
            ch = changes[0]
            ctype = ch.get("type")
            start = int(ch.get("start", 1))
            if ctype == "insert":
                end = start - 1
                inserted = str(ch.get("sequence", ""))
            elif ctype == "delete":
                end = int(ch.get("end", start))
                inserted = ""
            elif ctype == "delete_insert":
                start = start + 1
                end = int(ch.get("end", start)) - 1
                inserted = str(ch.get("sequence", ""))
            else:
                end = start
                inserted = str(ch.get("sequence", ""))

    # 2. If not found in repeat_dict, check known canonical variants directly
    if start is None:
        if mut_name in KNOWN_VARIANTS:
            canon = mut_name
            citation = KNOWN_VARIANTS[mut_name]
            enriched["canonical_name"] = canon
            enriched["is_known_variant"] = True
            enriched["literature_citation"] = citation
            enriched["hgvs_cdna"] = format_hgvs_cdna(
                canon, "duplication" if "dup" in canon else "insertion"
            )
            enriched["repeat_relative_coordinate"] = format_repeat_relative_coordinate(
                repeat_idx, canon
            )
            enriched["transcript_coordinate"] = "transcript_coordinate_unresolved"
            enriched["confidence_tier"] = NomenclatureRecord.TIER_A
            return enriched

        if "dup" in mut_name:
            event = "duplication"
        elif "delins" in mut_name:
            event = "delins"
        elif "del" in mut_name:
            event = "deletion"
        elif "ins" in mut_name:
            event = "insertion"
        else:
            event = "substitution"
        enriched["canonical_name"] = mut_name
        enriched["event_type"] = event
        enriched["hgvs_cdna"] = format_hgvs_cdna(mut_name, event)
        enriched["repeat_relative_coordinate"] = format_repeat_relative_coordinate(
            repeat_idx, mut_name
        )
        enriched["transcript_coordinate"] = "transcript_coordinate_unresolved"
        enriched["confidence_tier"] = (
            NomenclatureRecord.TIER_B if enriched.get("frameshift") else NomenclatureRecord.TIER_C
        )
        enriched["is_known_variant"] = False
        enriched["literature_citation"] = None
        enriched["repeat_form"] = None
        enriched["ambiguity_interval"] = None
        return enriched

    assert start is not None and end is not None and inserted is not None
    rec = name_variant_call(
        start,
        end,
        inserted,
        unit_symbol=closest_type,
        unit_sequence=CANONICAL_UNIT,
        support_reads=support_reads,
        independent_sources=1,
        repeat_index=repeat_idx,
    )
    rec_dict = rec.to_dict()
    enriched.update(rec_dict)
    return enriched
