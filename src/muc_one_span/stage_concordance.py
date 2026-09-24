"""Clair3 pileup-stage versus applied-call frameshift concordance.

Clair3 writes a pileup-stage VCF (``pileup.vcf.gz``) before its full-alignment
stage produces the final calls. When the pileup stage calls a frameshift indel
at a high allele fraction and depth that the applied (final) call set does not
contain, the absence of a frameshift in the consensus is not established.

This gate only prevents a reassuring negative. It never adds a variant to the
consensus and never turns a result positive; it records the discordant pileup
records so that callers can downgrade a NEGATIVE to INCONCLUSIVE or add a caveat.

Both thresholds come from ``CallingSettings``
(``stage_discordance_min_af`` and ``stage_discordance_min_depth``).

Normalisation: pileup records are split and left-aligned with
``bcftools norm -f <reference> -m -any``. Multi-ALT records in the final VCF are
split per genotype index but not re-trimmed, so an allele that is represented
differently can be reported as discordant. That errs towards over-flagging,
which is the safe direction for this gate.
"""

from __future__ import annotations

import contextlib
import math
import re
from pathlib import Path
from typing import Any

from muc_one_span.settings import CallingSettings
from muc_one_span.tools import run_tool
from muc_one_span.vcf import parse_vcf_variants

PILEUP_VCF_NAME = "pileup.vcf.gz"  # Clair3 output file name
STATUS_CONCORDANT = "concordant"
STATUS_DISCORDANT = "discordant_frameshift"
STATUS_NOT_ASSESSED = "not_assessed"
SOURCE = "clair3_pileup"
REASON_PILEUP_UNAVAILABLE = "pileup_vcf_unavailable"
REASON_FINAL_UNAVAILABLE = "final_vcf_unavailable"
NORMALIZED_PILEUP_NAME = "pileup.normalized.vcf.gz"
# bcftools query format: identity columns plus per-sample FORMAT AF and DP.
PILEUP_QUERY_FORMAT = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AF\\t%DP]\\n"
_QUERY_FIELDS = ("chrom", "pos", "ref", "alt", "af", "dp")
_SEQUENCE_ALLELE = re.compile(r"[ACGTNacgtn]+")
_GENOTYPE_SEPARATOR = re.compile(r"[/|]")

Allele = tuple[str, int, str, str]


def is_frameshift(ref: str, alt: str) -> bool:
    """Return whether a sequence-resolved indel changes length by a non-multiple of 3.

    Symbolic, spanning-deletion or otherwise non-sequence alleles are never frameshifts.
    """
    if not _SEQUENCE_ALLELE.fullmatch(ref) or not _SEQUENCE_ALLELE.fullmatch(alt):
        return False
    return (len(alt) - len(ref)) % 3 != 0


def applied_alleles(variants: list[dict]) -> set[Allele]:
    """Return ``(chrom, pos, ref, alt_i)`` for every nonzero genotype index in each record."""
    applied: set[Allele] = set()
    for variant in variants:
        alts = str(variant["alt"]).split(",")
        for index in _GENOTYPE_SEPARATOR.split(str(variant.get("genotype", "."))):
            if index not in {".", "0"}:
                alt = alts[int(index) - 1]
                applied.add((variant["chrom"], int(variant["pos"]), variant["ref"], alt))
    return applied


def discordant_frameshifts(
    pileup: list[dict], applied: set[Allele], *, min_af: float, min_depth: int
) -> list[dict]:
    """Return pileup frameshifts at ``af >= min_af`` and ``dp >= min_depth`` not applied."""
    records = [
        row
        for row in pileup
        if is_frameshift(row["ref"], row["alt"])
        and row["af"] >= min_af
        and row["dp"] >= min_depth
        and (row["chrom"], row["pos"], row["ref"], row["alt"]) not in applied
    ]
    return sorted(records, key=lambda row: (row["chrom"], row["pos"]))


def assess(pileup: list[dict], applied: set[Allele], settings: CallingSettings) -> dict[str, Any]:
    """Build the stage-concordance record for one allele's Clair3 partition."""
    records = discordant_frameshifts(
        pileup,
        applied,
        min_af=settings.stage_discordance_min_af,
        min_depth=settings.stage_discordance_min_depth,
    )
    return {
        "status": STATUS_DISCORDANT if records else STATUS_CONCORDANT,
        "source": SOURCE,
        "min_af": settings.stage_discordance_min_af,
        "min_depth": settings.stage_discordance_min_depth,
        "records": records,
    }


def _parse_float(value: str) -> float | None:
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _parse_int(value: str) -> int | None:
    try:
        return int(value)
    except ValueError:
        return None


def read_pileup(pileup_vcf: Path, reference: Path, work_dir: Path) -> list[dict]:
    """Normalise and query a Clair3 pileup VCF into one row per ALT allele.

    Rows whose AF or DP is missing or non-numeric are skipped because they cannot
    establish anything. Malformed query lines raise ``ValueError``.
    """
    normalized = work_dir / NORMALIZED_PILEUP_NAME
    try:
        run_tool(
            [
                "bcftools",
                "norm",
                "-f",
                str(reference),
                "-m",
                "-any",
                "-O",
                "z",
                "-o",
                str(normalized),
                str(pileup_vcf),
            ]
        )
        output = run_tool(["bcftools", "query", "-f", PILEUP_QUERY_FORMAT, str(normalized)])
    finally:
        with contextlib.suppress(FileNotFoundError):
            normalized.unlink()
    rows: list[dict] = []
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) != len(_QUERY_FIELDS):
            raise ValueError(f"Malformed pileup query record: {line!r}")
        chrom, pos, ref, alt, af_text, dp_text = fields
        af, dp = _parse_float(af_text), _parse_int(dp_text)
        if af is None or dp is None:
            continue
        rows.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, "af": af, "dp": dp})
    return rows


def _not_assessed(reason: str, settings: CallingSettings) -> dict[str, Any]:
    return {
        "status": STATUS_NOT_ASSESSED,
        "reason": reason,
        "source": SOURCE,
        "min_af": settings.stage_discordance_min_af,
        "min_depth": settings.stage_discordance_min_depth,
        "records": [],
    }


def annotate_stage_concordance(
    pileup_vcf: Path,
    final_vcf: Path,
    reference: Path,
    settings: CallingSettings | None = None,
    *,
    work_dir: Path | None = None,
) -> dict[str, Any]:
    """Compare a Clair3 pileup VCF with the applied final VCF of the same partition.

    Missing inputs return ``not_assessed`` without running any tool. ``work_dir``
    defaults to the pileup VCF's directory.
    """
    settings = settings or CallingSettings()
    if not pileup_vcf.is_file():
        return _not_assessed(REASON_PILEUP_UNAVAILABLE, settings)
    if not final_vcf.is_file():
        return _not_assessed(REASON_FINAL_UNAVAILABLE, settings)
    pileup = read_pileup(pileup_vcf, reference, work_dir or pileup_vcf.parent)
    return assess(pileup, applied_alleles(parse_vcf_variants(final_vcf)), settings)


def stage_concordance_reasons(info: Any, label: str) -> list[str]:
    """Return the report reason for a discordant allele, or ``[]`` otherwise."""
    record = info.get("stage_concordance") if isinstance(info, dict) else None
    if not isinstance(record, dict) or record.get("status") != STATUS_DISCORDANT:
        return []
    records = record["records"]
    first = records[0]
    return [
        f"{label}: caller-stage discordance: {len(records)} frameshift indel(s) called by "
        f"the Clair3 pileup stage at allele fraction >= {record['min_af']} and depth >= "
        f"{record['min_depth']} are not in the applied calls (first: {first['chrom']}:"
        f"{first['pos']} {first['ref']}>{first['alt']}, AF {first['af']}, DP {first['dp']}); "
        "absence of a frameshift is not established."
    ]
