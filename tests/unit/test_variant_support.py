"""Variant concordance requires exact event identity and replayed coordinates."""

import json
from pathlib import Path

import pytest

from muc_one_span.classify import classify_sequence, validate_mutations_against_vcf
from muc_one_span.config import load_repeat_dictionary


def fixture(tmp_path: Path, qual: float = 25) -> tuple:
    rd = load_repeat_dictionary()
    x = rd.repeats["X"]
    seq = x + x[:40] + "TT" + x[40:] + x
    ref = "ACGT" + x * 3 + "TGCA"
    full = "ACGT" + seq + "TGCA"
    rp, cp = tmp_path / "ref.fa", tmp_path / "consensus.fa"
    rp.write_text(">contig_1\n" + ref + "\n")
    cp.write_text(">contig_1\n" + full + "\n")
    pos = 4 + 60 + 40
    variants = [
        {
            "chrom": "contig_1",
            "pos": pos,
            "ref": ref[pos - 1],
            "alt": ref[pos - 1] + "TT",
            "genotype": "1/1",
            "qual": qual,
        }
    ]
    context = {
        "reference_path": str(rp),
        "full_consensus_path": str(cp),
        "trim_start": 4,
        "trim_end": len(full) - 4,
        "chrom": "contig_1",
        "haplotype": "I",
        "sample": "S",
    }
    return rd, seq, variants, context


def test_exact_indel_concordance(tmp_path: Path) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    result = validate_mutations_against_vcf(
        classify_sequence(seq, rd),
        variants,
        sequence=seq,
        repeat_dict=rd,
        consensus_context=context,
    )
    assert result["mutations_detected"][0]["vcf_support"] is True
    assert result["mutations_detected"][0]["vcf_support_status"] == "exact_sequence_concordance"


@pytest.mark.parametrize(
    "change",
    [{"alt": "C"}, {"chrom": "wrong"}, {"genotype": "0/0"}, {"genotype": "0/1"}, {"alt": "<INS>"}],
)
def test_other_or_unselected_event_cannot_support(tmp_path: Path, change: dict) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    variants[0].update(change)
    result = validate_mutations_against_vcf(
        classify_sequence(seq, rd),
        variants,
        sequence=seq,
        repeat_dict=rd,
        consensus_context=context,
    )
    assert result["mutations_detected"][0]["vcf_support"] is False


def test_upstream_flank_indel_projects_to_consensus(tmp_path: Path) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    cp = Path(context["full_consensus_path"])
    full = cp.read_text().splitlines()[1]
    cp.write_text(">contig_1\n" + full[:1] + "AA" + full[1:] + "\n")
    context["trim_start"] += 2
    context["trim_end"] += 2
    variants.insert(
        0, {"chrom": "contig_1", "pos": 1, "ref": "A", "alt": "AAA", "genotype": "1/1", "qual": 30}
    )
    result = validate_mutations_against_vcf(
        classify_sequence(seq, rd),
        variants,
        sequence=seq,
        repeat_dict=rd,
        consensus_context=context,
    )
    assert result["mutations_detected"][0]["vcf_support"] is True
    # Context is a JSON-serializable artifact contract.
    json.dumps(context)


def test_wrong_trim_disables_projection(tmp_path: Path) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    context["trim_start"] += 1
    result = validate_mutations_against_vcf(
        classify_sequence(seq, rd),
        variants,
        sequence=seq,
        repeat_dict=rd,
        consensus_context=context,
    )
    assert result["mutations_detected"][0]["vcf_support"] is False


@pytest.mark.parametrize(
    "problem,reason",
    [
        ("multiallelic_heterozygous_indel", "ambiguous_genotype_selection"),
        ("trim", "trim_mismatch"),
        ("consensus", "consensus_replay_mismatch"),
    ],
)
def test_projection_exposes_allele_failure_reason(tmp_path, problem, reason):
    rd, seq, variants, context = fixture(tmp_path)
    if problem == "multiallelic_heterozygous_indel":
        anchor = variants[0]["ref"]
        variants[0].update(alt=f"{anchor}TT,{anchor}T", genotype="1/2")
    elif problem == "trim":
        context["trim_start"] += 1
    else:
        path = Path(context["full_consensus_path"])
        path.write_text(path.read_text() + "A\n")
    result = validate_mutations_against_vcf(
        classify_sequence(seq, rd),
        variants,
        sequence=seq,
        repeat_dict=rd,
        consensus_context=context,
    )
    assert result["vcf_projection"]["status"] == "unavailable"
    assert result["vcf_projection"]["reason"] == reason


def _validate(rd, seq: str, variants: list[dict], context: dict) -> dict:
    return validate_mutations_against_vcf(
        classify_sequence(seq, rd),
        variants,
        sequence=seq,
        repeat_dict=rd,
        consensus_context=context,
    )


def test_unrelated_heterozygous_indel_keeps_per_event_support(tmp_path: Path) -> None:
    """bcftools -H I applies a 0/1 flank insertion; the homozygous VNTR event stays supported."""
    rd, seq, variants, context = fixture(tmp_path)
    reference = Path(context["reference_path"]).read_text().splitlines()[1]
    flank = len(reference) - 4  # right flank "TGCA"; insert A after its G
    assert reference[flank + 1] == "G"
    variants.append(
        {"chrom": "contig_1", "pos": flank + 2, "ref": "G", "alt": "GA", "genotype": "0/1"}
    )
    consensus = Path(context["full_consensus_path"])
    full = consensus.read_text().splitlines()[1]
    consensus.write_text(">contig_1\n" + full[:-2] + "A" + full[-2:] + "\n")
    context["trim_end"] = len(full) + 1 - 5
    result = _validate(rd, seq, variants, context)
    mutation = result["mutations_detected"][0]
    assert result["vcf_projection"]["status"] == "available"
    assert result["vcf_projection"]["unresolved_genotype_edits"] == 1
    assert mutation["vcf_support_status"] == "exact_sequence_concordance"
    assert mutation["vcf_support"] is True


def test_heterozygous_event_under_iupac_is_unresolved_not_supported(tmp_path: Path) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    variants[0]["genotype"] = "0/1"
    result = _validate(rd, seq, variants, context)
    mutation = result["mutations_detected"][0]
    assert result["vcf_projection"]["status"] == "available"
    assert mutation["vcf_support_status"] == "heterozygous_genotype_unresolved"
    assert mutation["vcf_support"] is False


def test_phase_selected_heterozygous_event_remains_supported(tmp_path: Path) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    variants[0]["genotype"] = "0|1"
    context["haplotype"] = 2
    mutation = _validate(rd, seq, variants, context)["mutations_detected"][0]
    assert mutation["vcf_support_status"] == "exact_sequence_concordance"
