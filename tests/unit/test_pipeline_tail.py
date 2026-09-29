"""``finish_run``'s classification-confidence echo for the hybrid engine.

The ladder's ``classify.py`` ``allele_confidence`` is a dictionary-fit heuristic fed
identically for both engines (see ``hybrid/allele_fields.py::allele_info``), so it
carries no hybrid reconstruction evidence. ``finish_run`` additively echoes the hybrid
engine's own read-support evidence (``consensus_concordance_fraction``, present only
on hybrid allele records) alongside the existing ``confidence:`` line, never in its
place.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.pipeline_tail import finish_run
from muc_one_span.settings import DEFAULT_SETTINGS

RD = load_repeat_dictionary()
FAKE_CLASSIFY = {
    "structure": "1-2-3-X-6-7-8-9",
    "repeats": [],
    "mutations_detected": [],
    "allele_confidence": 1.0,
}


def _finish(tmp_path: Path, alleles_result: dict) -> str:
    fa = tmp_path / "allele_1.fa"
    fa.write_text(">allele_1\nACGT\n")
    with (
        patch("muc_one_span.classify.classify_sequence", return_value=FAKE_CLASSIFY),
        patch("muc_one_span.classify.validate_mutations_against_vcf", return_value=FAKE_CLASSIFY),
    ):
        summary = finish_run(
            out=tmp_path,
            input_path="in.fastq",
            rd=RD,
            settings=replace(DEFAULT_SETTINGS, run=replace(DEFAULT_SETTINGS.run, engine="hybrid")),
            alleles_result=alleles_result,
            consensus_paths={"allele_1": fa},
            vcf_paths={},
            tool_versions={},
            configuration_record={},
            report=False,
            bam_path=None,
            fasta_path=None,
        )
    assert summary["run_status"]["status"] == "completed"
    return summary


def test_hybrid_agreement_is_echoed_alongside_the_dictionary_confidence(
    tmp_path: Path, capsys
) -> None:
    _finish(tmp_path, {"allele_1": {"consensus_concordance_fraction": 0.42}})
    out = capsys.readouterr().out
    assert "confidence: 1.00" in out
    assert "consensus concordance (hybrid read support): 0.42" in out


def test_ladder_allele_without_the_hybrid_field_gets_no_extra_line(tmp_path: Path, capsys) -> None:
    _finish(tmp_path, {"allele_1": {"length": 30}})
    out = capsys.readouterr().out
    assert "confidence: 1.00" in out
    assert "consensus concordance" not in out
