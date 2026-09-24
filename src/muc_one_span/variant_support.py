"""Exact sequence concordance with the VCF used to construct a consensus.

This is internal consistency evidence, not independent validation. Projection is
available only when selected alleles replay the complete consensus exactly.
"""

from __future__ import annotations

from pathlib import Path

from muc_one_span.config import RepeatDictionary

_IUPAC = {
    frozenset("AC"): "M",
    frozenset("AG"): "R",
    frozenset("AT"): "W",
    frozenset("CG"): "S",
    frozenset("CT"): "Y",
    frozenset("GT"): "K",
    frozenset("ACG"): "V",
    frozenset("ACT"): "H",
    frozenset("AGT"): "D",
    frozenset("CGT"): "B",
    frozenset("ACGT"): "N",
}


def _fasta(path: str) -> tuple[str, str]:
    lines = Path(path).read_text().splitlines()
    headers = [line[1:].split()[0] for line in lines if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError("Concordance context requires one FASTA record")
    return headers[0], "".join(line.strip() for line in lines if not line.startswith(">"))


def _selected_alt(variant: dict, haplotype: str | int) -> str | None:
    alleles: list[str] = [variant["ref"], *variant["alt"].split(",")]
    gt = variant["genotype"].replace("|", "/").split("/")
    if len(gt) != 2 or any(not g.isdigit() or int(g) >= len(alleles) for g in gt):
        return None
    selected = [alleles[int(g)] for g in gt]
    if haplotype in (1, 2, "1", "2"):
        return selected[int(haplotype) - 1]
    if haplotype != "I":
        return None
    if selected[0] == selected[1]:
        return selected[0]
    # SNP IUPAC has a precise coordinate map; indel ambiguity does not.
    if len(variant["ref"]) == 1 and all(len(a) == 1 and a in "ACGT" for a in selected):
        return _IUPAC.get(frozenset(selected))
    # bcftools consensus -H I applies the ALT of a heterozygous REF/ALT indel
    # (verified with bcftools 1.17; tests/integration guards the installed tool).
    # Multi-ALT heterozygous indels and MNPs have no verified mapping.
    alts = [a for a in selected if a != variant["ref"]]
    if len(alts) == 1 and len(alts[0]) != len(variant["ref"]):
        return alts[0]
    return None


def _replay(variants: list[dict], context: dict) -> tuple[str, list[dict]] | str:
    chrom, reference = _fasta(context["reference_path"])
    consensus_chrom, consensus = _fasta(context["full_consensus_path"])
    if chrom != context["chrom"] or consensus_chrom != chrom:
        return "contig_mismatch"
    pieces: list[str] = []
    applied: list[dict] = []
    ref_pos, query_pos = 0, 0
    for variant in sorted(variants, key=lambda v: v["pos"]):
        if variant.get("chrom") != chrom:
            return "variant_contig_mismatch"
        ref = variant.get("ref", "")
        if not ref or any(b not in "ACGT" for b in ref):
            return "unsupported_reference_allele"
        alt = _selected_alt(variant, context.get("haplotype", "I"))
        if alt is None:
            return "ambiguous_genotype_selection"
        if not alt or any(b not in "ACGTMRWSYKVHDBN" for b in alt):
            return "unsupported_alternate_allele"
        start = variant["pos"] - 1
        if start < ref_pos:
            return "overlapping_variant_records"
        if reference[start : start + len(ref)] != ref:
            return "reference_allele_mismatch"
        unchanged = reference[ref_pos:start]
        pieces.extend((unchanged, alt))
        query_pos += len(unchanged)
        genotype = set(variant["genotype"].replace("|", "/").split("/"))
        unresolved = context.get("haplotype", "I") == "I" and len(genotype) > 1
        applied.append(
            {
                **variant,
                "selected_alt": alt,
                "query_start": query_pos,
                "unresolved_genotype": unresolved,
            }
        )
        query_pos += len(alt)
        ref_pos = start + len(ref)
    pieces.append(reference[ref_pos:])
    return (consensus, applied) if "".join(pieces) == consensus else "consensus_replay_mismatch"


def mutation_concordance(
    classification: dict,
    variants: list[dict],
    sequence: str | None,
    repeat_dict: RepeatDictionary | None,
    context: dict | None,
) -> dict[int, tuple[str, list[dict]]]:
    """Match a repeat event by reverting its selected VCF allele in full context.

    Whole-sequence comparison permits equivalent homopolymer anchoring without
    confusing a nearby unrelated SNP/indel with the classified event. If reverting
    the same VCF event explains multiple reported repeats, localization is ambiguous.
    Unsupported complex/missing context is explicit, never proximity support.
    """
    mutations = classification.get("mutations_detected", [])
    result: dict[int, tuple[str, list[dict]]] = {
        m["repeat_index"]: ("projection_unavailable", []) for m in mutations
    }
    projection = {"status": "unavailable", "reason": "missing_context"}
    classification["vcf_projection"] = projection
    if context is None or sequence is None or repeat_dict is None:
        return result
    required = {"chrom", "pos", "ref", "alt", "genotype"}
    if any(not required.issubset(v) for v in variants):
        projection["reason"] = "incomplete_variant_fields"
        return result
    replayed = _replay(variants, context)
    if isinstance(replayed, str):
        projection["reason"] = replayed
        return result
    full, edits = replayed
    trim_start, trim_end = context["trim_start"], context["trim_end"]
    if not 0 <= trim_start <= trim_end <= len(full) or full[trim_start:trim_end] != sequence:
        projection["reason"] = "trim_mismatch"
        return result
    projection.update(
        status="available",
        reason="exact_full_consensus_replay",
        unresolved_genotype_edits=sum(edit["unresolved_genotype"] for edit in edits),
    )
    used: dict[int, list[int]] = {}
    for mutation in mutations:
        index = mutation["repeat_index"]
        repeat = classification["repeats"][index - 1]
        if repeat.get("localization_status", "resolved") != "resolved":
            result[index] = ("localization_ambiguous", [])
            continue
        parent = repeat_dict.repeats.get(mutation["closest_type"])
        if parent is None or "start" not in repeat or "end" not in repeat:
            continue
        start, end = trim_start + repeat["start"], trim_start + repeat["end"]
        reverted_repeat = full[:start] + parent + full[end:]
        supporting: list[dict] = []
        for edit_index, edit in enumerate(edits):
            alt, ref, pos = edit["selected_alt"], edit["ref"], edit["query_start"]
            if len(alt) == len(ref) or any(b not in "ACGT" for b in alt):
                continue
            reverted_event = full[:pos] + ref + full[pos + len(alt) :]
            if reverted_event == reverted_repeat:
                supporting.append(edit)
                used.setdefault(edit_index, []).append(index)
        resolved = [edit for edit in supporting if not edit["unresolved_genotype"]]
        if resolved:
            result[index] = ("exact_sequence_concordance", resolved)
        elif supporting:
            result[index] = ("heterozygous_genotype_unresolved", [])
        else:
            result[index] = ("absent", [])
    for indices in used.values():
        if len(set(indices)) > 1:
            for index in indices:
                result[index] = ("localization_ambiguous", [])
    return result
