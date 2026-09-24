"""Clair3 pileup-stage versus applied-call frameshift concordance (issue #72)."""

from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span import stage_concordance as sc
from muc_one_span.settings import CallingSettings

R9_ROW = {"chrom": "contig_82", "pos": 2833, "ref": "C", "alt": "CG", "af": 0.719, "dp": 121}
DEFAULTS = CallingSettings()


def _row(**overrides: Any) -> dict[str, Any]:
    return {**R9_ROW, **overrides}


def _discordant(rows: list[dict[str, Any]], applied: set | None = None) -> list[dict]:
    return sc.discordant_frameshifts(
        rows,
        applied or set(),
        min_af=DEFAULTS.stage_discordance_min_af,
        min_depth=DEFAULTS.stage_discordance_min_depth,
    )


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("C", "CG", True),
        ("CCCCAG", "C", True),
        ("C", "CGGG", False),
        ("GCC", "G", True),
        ("C", "G", False),
        ("C", "*", False),
        ("C", "<DEL>", False),
    ],
)
def test_is_frameshift(ref: str, alt: str, expected: bool) -> None:
    assert sc.is_frameshift(ref, alt) is expected


def _variant(genotype: str, alt: str = "CG,CGG") -> dict[str, Any]:
    return {"chrom": "contig_82", "pos": 2833, "ref": "C", "alt": alt, "genotype": genotype}


@pytest.mark.parametrize(
    ("genotype", "expected"),
    [
        ("1/1", {("contig_82", 2833, "C", "CG")}),
        ("1/2", {("contig_82", 2833, "C", "CG"), ("contig_82", 2833, "C", "CGG")}),
        ("0/0", set()),
        ("./.", set()),
        (".", set()),
        ("0|1", {("contig_82", 2833, "C", "CG")}),
    ],
)
def test_applied_alleles(genotype: str, expected: set) -> None:
    assert sc.applied_alleles([_variant(genotype)]) == expected


def test_unapplied_r9_frameshift_is_discordant() -> None:
    assert _discordant([R9_ROW]) == [R9_ROW]


def test_applied_frameshift_is_not_discordant() -> None:
    assert _discordant([R9_ROW], {("contig_82", 2833, "C", "CG")}) == []


@pytest.mark.parametrize(
    ("overrides", "included"),
    [
        ({"af": 0.5}, True),
        ({"af": 0.4999}, False),
        ({"dp": 10}, True),
        ({"dp": 9}, False),
        ({"af": 0.156, "dp": 128}, False),
        ({"alt": "CGGG", "af": 0.9}, False),
        ({"alt": "A", "af": 0.745}, False),
    ],
)
def test_discordance_thresholds_and_shapes(overrides: dict[str, Any], included: bool) -> None:
    row = _row(**overrides)
    assert _discordant([row]) == ([row] if included else [])


def test_discordant_records_are_sorted_by_position() -> None:
    later, earlier = _row(pos=3000), _row(pos=100)
    assert _discordant([later, earlier]) == [earlier, later]


def test_assess_echoes_settings_and_status() -> None:
    settings = CallingSettings(stage_discordance_min_af=0.6, stage_discordance_min_depth=20)
    record = sc.assess([R9_ROW], set(), settings)
    assert record == {
        "status": sc.STATUS_DISCORDANT,
        "source": "clair3_pileup",
        "min_af": 0.6,
        "min_depth": 20,
        "records": [R9_ROW],
    }
    concordant = sc.assess([R9_ROW], {("contig_82", 2833, "C", "CG")}, settings)
    assert concordant["status"] == sc.STATUS_CONCORDANT
    assert concordant["records"] == []


def test_assess_honours_nondefault_thresholds() -> None:
    strict = CallingSettings(stage_discordance_min_af=0.8)
    assert sc.assess([R9_ROW], set(), strict)["status"] == sc.STATUS_CONCORDANT


def test_missing_pileup_is_not_assessed_without_tools(tmp_path: Path) -> None:
    final = tmp_path / "final.vcf.gz"
    final.touch()
    with patch("muc_one_span.stage_concordance.run_tool") as run_tool:
        record = sc.annotate_stage_concordance(
            tmp_path / sc.PILEUP_VCF_NAME, final, tmp_path / "ref.fa"
        )
    run_tool.assert_not_called()
    assert record["status"] == sc.STATUS_NOT_ASSESSED
    assert record["reason"] == "pileup_vcf_unavailable"


def test_missing_final_vcf_is_not_assessed_without_tools(tmp_path: Path) -> None:
    pileup = tmp_path / sc.PILEUP_VCF_NAME
    pileup.touch()
    with (
        patch("muc_one_span.stage_concordance.run_tool") as run_tool,
        patch("muc_one_span.stage_concordance.parse_vcf_variants") as parse,
    ):
        record = sc.annotate_stage_concordance(
            pileup, tmp_path / "final.vcf.gz", tmp_path / "ref.fa"
        )
    run_tool.assert_not_called()
    parse.assert_not_called()
    assert record["status"] == sc.STATUS_NOT_ASSESSED
    assert record["reason"] == "final_vcf_unavailable"


def test_read_pileup_normalizes_queries_and_cleans_up(tmp_path: Path) -> None:
    pileup, reference = tmp_path / sc.PILEUP_VCF_NAME, tmp_path / "ref.fa"
    normalized = tmp_path / "pileup.normalized.vcf.gz"
    query_output = "contig_82\t2833\tC\tCG\t0.719\t121\ncontig_82\t10\tC\tA\t.\t5\n"

    def fake_run(cmd: list[str], **_: Any) -> str:
        if cmd[1] == "norm":
            Path(cmd[cmd.index("-o") + 1]).write_text("normalized")
            return ""
        assert normalized.exists()
        return query_output

    with patch("muc_one_span.stage_concordance.run_tool", side_effect=fake_run) as run_tool:
        rows = sc.read_pileup(pileup, reference, tmp_path)
    first, second = (call.args[0] for call in run_tool.call_args_list)
    assert first[:6] == ["bcftools", "norm", "-f", str(reference), "-m", "-any"]
    assert first[-3:] == ["-o", str(normalized), str(pileup)]
    assert second[:2] == ["bcftools", "query"]
    assert second[-1] == str(normalized)
    assert rows == [R9_ROW]
    assert not normalized.exists()


def test_read_pileup_cleans_up_when_query_fails(tmp_path: Path) -> None:
    normalized = tmp_path / "pileup.normalized.vcf.gz"

    def fake_run(cmd: list[str], **_: Any) -> str:
        if cmd[1] == "norm":
            normalized.write_text("normalized")
            return ""
        raise RuntimeError("query failed")

    with (
        patch("muc_one_span.stage_concordance.run_tool", side_effect=fake_run),
        pytest.raises(RuntimeError, match="query failed"),
    ):
        sc.read_pileup(tmp_path / sc.PILEUP_VCF_NAME, tmp_path / "ref.fa", tmp_path)
    assert not normalized.exists()


def test_read_pileup_skips_nonnumeric_depth_and_rejects_malformed(tmp_path: Path) -> None:
    with patch(
        "muc_one_span.stage_concordance.run_tool", side_effect=["", "c\t5\tC\tCG\t0.9\tx\n"]
    ):
        assert sc.read_pileup(tmp_path / "p.vcf.gz", tmp_path / "r.fa", tmp_path) == []
    with (
        patch("muc_one_span.stage_concordance.run_tool", side_effect=["", "c\t5\tC\n"]),
        pytest.raises(ValueError, match="Malformed"),
    ):
        sc.read_pileup(tmp_path / "p.vcf.gz", tmp_path / "r.fa", tmp_path)


def test_annotate_end_to_end_is_discordant(tmp_path: Path) -> None:
    pileup, final = tmp_path / sc.PILEUP_VCF_NAME, tmp_path / "final.vcf.gz"
    pileup.touch()
    final.touch()
    reference = tmp_path / "ref.fa"
    with (
        patch("muc_one_span.stage_concordance.read_pileup", return_value=[R9_ROW]) as read,
        patch("muc_one_span.stage_concordance.parse_vcf_variants", return_value=[]),
    ):
        record = sc.annotate_stage_concordance(pileup, final, reference)
    read.assert_called_once_with(pileup, reference, pileup.parent)
    assert record["status"] == sc.STATUS_DISCORDANT
    assert record["records"] == [R9_ROW]
    assert (record["min_af"], record["min_depth"]) == (0.5, 10)


@pytest.mark.parametrize(
    "info",
    [
        {"stage_concordance": {"status": "concordant", "records": []}},
        {"stage_concordance": {"status": "not_assessed", "reason": "pileup_vcf_unavailable"}},
        {},
        None,
        "discordant_frameshift",
    ],
)
def test_stage_concordance_reasons_silent_unless_discordant(info: Any) -> None:
    assert sc.stage_concordance_reasons(info, "allele_1") == []


def test_stage_concordance_reason_text_is_exact() -> None:
    info = {"stage_concordance": sc.assess([R9_ROW], set(), DEFAULTS)}
    assert sc.stage_concordance_reasons(info, "allele_1") == [
        "allele_1: caller-stage discordance: 1 frameshift indel(s) called by the Clair3 "
        "pileup stage at allele fraction >= 0.5 and depth >= 10 are not in the applied "
        "calls (first: contig_82:2833 C>CG, AF 0.719, DP 121); absence of a frameshift "
        "is not established."
    ]


def test_discordant_status_without_records_still_blocks() -> None:
    """A hand-edited discordant record with no rows must not raise and must still block."""
    info = {"stage_concordance": {"status": "discordant_frameshift", "records": []}}
    reasons = sc.stage_concordance_reasons(info, "Allele 1")
    assert len(reasons) == 1
    assert reasons[0].startswith("Allele 1: caller-stage discordance:")
    assert reasons[0].endswith("absence of a frameshift is not established.")


@pytest.mark.parametrize("reason", ["pileup_vcf_unavailable", "final_vcf_unavailable"])
def test_not_assessed_yields_quality_caveat(reason: str) -> None:
    info = {"stage_concordance": {"status": "not_assessed", "reason": reason}}
    assert sc.stage_concordance_caveats(info, "Allele 2") == [
        f"Allele 2: caller-stage concordance not assessed ({reason}); Clair3 pileup-stage "
        "frameshift calls were not compared with the applied calls."
    ]


@pytest.mark.parametrize(
    "info",
    [
        {"stage_concordance": {"status": "concordant", "records": []}},
        {"stage_concordance": {"status": "discordant_frameshift", "records": [R9_ROW]}},
        {},
        None,
    ],
)
def test_stage_concordance_caveats_only_for_not_assessed(info: Any) -> None:
    assert sc.stage_concordance_caveats(info, "Allele 1") == []


def test_missing_pileup_logs_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    record = sc.annotate_stage_concordance(
        tmp_path / sc.PILEUP_VCF_NAME, tmp_path / "final.vcf.gz", tmp_path / "ref.fa"
    )
    assert record["status"] == sc.STATUS_NOT_ASSESSED
    assert "pileup_vcf_unavailable" in caplog.text
    assert caplog.records[-1].levelname == "WARNING"
