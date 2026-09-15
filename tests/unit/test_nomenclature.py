"""Unit tests for MUC1 VNTR HGVS nomenclature and 3'-most normalization."""

from __future__ import annotations

from muc_one_span.nomenclature import (
    CANONICAL_UNIT,
    NomenclatureRecord,
    ambiguity_interval,
    enrich_mutation_record,
    name_edit,
    name_variant_call,
    normalise,
    repeat_form,
    revcomp,
    trim_edit,
)


def test_revcomp() -> None:
    """Verify DNA reverse complementation including ambiguous N."""
    assert revcomp("ACGTN") == "NACGT"
    assert revcomp("GCCCACGGTG") == "CACCGTGGGC"


def test_trim_edit_shared_flanks() -> None:
    """Verify shared prefix and suffix trimming."""
    # Prefix match
    s, e, ins = trim_edit(CANONICAL_UNIT, 1, 2, "GA")
    assert (s, e, ins) == (2, 2, "A")

    # Suffix match
    s, e, ins = trim_edit(CANONICAL_UNIT, 1, 2, "TC")
    assert (s, e, ins) == (1, 1, "T")


def test_canonical_59dupc_normalization() -> None:
    """An insertion of C anywhere in the 7xC tract (53-59) normalizes to 59dupC."""
    # Insertion at 53 (54, 53)
    canon_name, event = name_edit(CANONICAL_UNIT, 54, 53, "C")
    assert canon_name == "59dupC"
    assert event == "duplication"

    # Insertion at 56 (57, 56)
    canon_name, event = name_edit(CANONICAL_UNIT, 57, 56, "C")
    assert canon_name == "59dupC"
    assert event == "duplication"

    # Insertion at 59 (60, 59)
    canon_name, event = name_edit(CANONICAL_UNIT, 60, 59, "C")
    assert canon_name == "59dupC"
    assert event == "duplication"


def test_canonical_59dupc_ambiguity_interval() -> None:
    """The 7xC homopolymer tract yields an ambiguity interval of (53, 59)."""
    interval = ambiguity_interval(CANONICAL_UNIT, 54, 53, "C")
    assert interval == (53, 59)


def test_canonical_59dupc_repeat_form() -> None:
    """A single C duplication in 7xC tract converts to 53C[7]>53C[8]."""
    rep = repeat_form(CANONICAL_UNIT, 54, 53, "C")
    assert rep == "53C[7]>53C[8]"


def test_multi_base_56_59dupcccc() -> None:
    """An insertion of 4 Cs normalizes to 56_59dupCCCC and 53C[7]>53C[11]."""
    canon_name, event = name_edit(CANONICAL_UNIT, 54, 53, "CCCC")
    assert canon_name == "56_59dupCCCC"
    assert event == "duplication"
    assert repeat_form(CANONICAL_UNIT, 54, 53, "CCCC") == "53C[7]>53C[11]"


def test_58_59insg_interruption() -> None:
    """An inserted G after position 58 cannot be dup; names as 58_59insG."""
    # Reference bases 58=C, 59=C, 60=A
    canon_name, event = name_edit(CANONICAL_UNIT, 59, 58, "G")
    assert canon_name == "58_59insG"
    assert event == "insertion"
    # Not a homopolymer run of Gs
    assert repeat_form(CANONICAL_UNIT, 59, 58, "G") is None


def test_60dupa_terminal_duplication() -> None:
    """Duplication of the terminal A at position 60 names as 60dupA."""
    canon_name, event = name_edit(CANONICAL_UNIT, 61, 60, "A")
    assert canon_name == "60dupA"
    assert event == "duplication"


def test_54_56delinsat_anchored() -> None:
    """Delins variants are anchored and never shifted."""
    canon_name, event = name_edit(CANONICAL_UNIT, 54, 56, "AT")
    assert canon_name == "54_56delinsAT"
    assert event == "delins"
    assert ambiguity_interval(CANONICAL_UNIT, 54, 56, "AT") is None


def test_deletion_1_5delgccca() -> None:
    """A 5 bp deletion at unit start names as 1_5delGCCCA."""
    canon_name, event = name_edit(CANONICAL_UNIT, 1, 5, "")
    assert canon_name == "1_5delGCCCA"
    assert event == "deletion"


def test_single_nucleotide_substitution() -> None:
    """A single nucleotide change names with standard > notation."""
    # Base 2 is C in GCCCACGGTG...
    canon_name, event = name_edit(CANONICAL_UNIT, 2, 2, "A")
    assert canon_name == "2C>A"
    assert event == "substitution"


def test_name_variant_call_record() -> None:
    """Full NomenclatureRecord generation with clinical fields."""
    rec = name_variant_call(54, 53, "C", support_reads=15, independent_sources=2)
    assert rec.canonical_name == "59dupC"
    assert rec.event_type == "duplication"
    assert rec.ambiguity_interval == (53, 59)
    assert rec.repeat_form == "53C[7]>53C[8]"
    assert rec.hgvs_cdna == "NM_001204286.1:c.59dupC"
    assert rec.confidence_tier == NomenclatureRecord.TIER_A
    assert rec.is_known_variant is True
    assert "Kirby et al. 2013" in str(rec.literature_citation)

    d = rec.to_dict()
    assert d["canonical_name"] == "59dupC"
    assert d["ambiguity_interval"] == [53, 59]


def test_low_support_tier_demotion() -> None:
    """Variants with thin read support are demoted to Tier B or Tier C."""
    # Support 3 is Tier B
    rec_b = name_variant_call(54, 53, "C", support_reads=3, independent_sources=1)
    assert rec_b.confidence_tier == NomenclatureRecord.TIER_B

    # Support 1 is Tier C
    rec_c = name_variant_call(54, 53, "C", support_reads=1, independent_sources=1)
    assert rec_c.confidence_tier == NomenclatureRecord.TIER_C


def test_deletion_rolling_and_ambiguity() -> None:
    """A single C deletion in the 7xC tract rolls 3' to 59delC with ambiguity (53, 59)."""
    canon_name, event = name_edit(CANONICAL_UNIT, 54, 54, "")
    assert canon_name == "59delC"
    assert event == "deletion"

    interval = ambiguity_interval(CANONICAL_UNIT, 54, 54, "")
    assert interval == (53, 59)

    rep = repeat_form(CANONICAL_UNIT, 54, 54, "")
    assert rep == "53C[7]>53C[6]"


def test_classify_event_and_formatting() -> None:
    """Test all event classification and HGVS formatting branches."""
    from muc_one_span.nomenclature import classify_event, format_hgvs_cdna, is_duplication

    # Duplication
    assert classify_event(60, 59, "C") == "duplication"
    # Insertion
    assert classify_event(59, 58, "G") == "insertion"
    # Deletion
    assert classify_event(10, 15, "") == "deletion"
    # Substitution
    assert classify_event(5, 5, "T") == "substitution"
    # Delins
    assert classify_event(5, 7, "AA") == "delins"

    # Formatting branches
    assert format_hgvs_cdna("59dupC", "duplication") == "NM_001204286.1:c.59dupC"
    assert format_hgvs_cdna("58_59insG", "insertion") == "NM_001204286.1:c.58_59insG"
    assert format_hgvs_cdna("1_5delGCCCA", "deletion") == "NM_001204286.1:c.1_5delGCCCA"
    assert format_hgvs_cdna("54_56delinsAT", "delins") == "NM_001204286.1:c.54_56delinsAT"
    assert format_hgvs_cdna("2C>A", "substitution") == "NM_001204286.1:c.2C>A"
    assert format_hgvs_cdna("other", "unknown") == "NM_001204286.1:c.other"

    # Duplication check when left < ins_len
    assert is_duplication(CANONICAL_UNIT, 1, "GCCC") is False


def test_repeat_form_edge_cases() -> None:
    """Test repeat form with mixed inserted sequence or mismatched base."""
    # Mixed insertion in tract
    assert repeat_form(CANONICAL_UNIT, 54, 53, "CA") is None
    # Mixed deletion
    assert repeat_form(CANONICAL_UNIT, 52, 54, "") is None


def test_out_of_bounds_handling() -> None:
    """Out-of-bounds start/end are returned without crash."""
    s, e, ins = normalise(CANONICAL_UNIT, 0, 10, "A")
    assert (s, e, ins) == (0, 10, "A")

    interval = ambiguity_interval(CANONICAL_UNIT, 0, 70, "A")
    assert interval is None

    rep = repeat_form(CANONICAL_UNIT, 0, 70, "A")
    assert rep is None


def test_enrich_mutation_record() -> None:
    """Enriching mutation dict resolves repeats.json definitions to HGVS nomenclature."""
    # Test known literature mutation dupC
    m_dupc = {
        "repeat_index": 45,
        "closest_type": "X",
        "mutation_name": "dupC",
        "frameshift": True,
        "vcf_support": True,
    }
    enriched = enrich_mutation_record(m_dupc)
    assert enriched["canonical_name"] == "59dupC"
    assert enriched["hgvs_cdna"] == "NM_001204286.1:c.59dupC"
    assert enriched["repeat_form"] == "53C[7]>53C[8]"
    assert enriched["ambiguity_interval"] == [53, 59]
    assert enriched["confidence_tier"] == NomenclatureRecord.TIER_A
    assert enriched["is_known_variant"] is True
    assert "Kirby et al." in str(enriched["literature_citation"])

    # Test known literature mutation insG
    m_insg = {
        "repeat_index": 20,
        "closest_type": "X",
        "mutation_name": "insG",
        "frameshift": True,
        "vcf_support": True,
    }
    enriched_insg = enrich_mutation_record(m_insg)
    assert enriched_insg["canonical_name"] == "58_59insG"
    assert enriched_insg["hgvs_cdna"] == "NM_001204286.1:c.58_59insG"
    assert enriched_insg["confidence_tier"] == NomenclatureRecord.TIER_A
    assert enriched_insg["is_known_variant"] is True

    # Test novel frameshift variant
    m_novel_fs = {
        "repeat_index": 12,
        "closest_type": "X",
        "mutation_name": "novel_ins1bp",
        "frameshift": True,
        "vcf_support": True,
    }
    enriched_fs = enrich_mutation_record(m_novel_fs)
    assert enriched_fs["confidence_tier"] == NomenclatureRecord.TIER_B
    assert enriched_fs["is_known_variant"] is False

    # Test novel in-frame variant
    m_novel_inframe = {
        "repeat_index": 12,
        "closest_type": "X",
        "mutation_name": "novel_sub",
        "frameshift": False,
        "vcf_support": False,
    }
    enriched_inframe = enrich_mutation_record(m_novel_inframe)
    assert enriched_inframe["confidence_tier"] == NomenclatureRecord.TIER_C
    assert enriched_inframe["is_known_variant"] is False
