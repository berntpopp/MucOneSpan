"""Deterministic synthetic MUC1 alleles and noisy long reads for hybrid-engine tests."""

from __future__ import annotations

import random

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.spans import ReadRecord

RD = load_repeat_dictionary()
PRE = ["1", "2", "3", "4", "5"]
POST = ["6", "7", "8", "9"]


def dupc(unit: str = "X") -> str:
    """Parent unit carrying the dictionary dupC template (7C->8C; insert at 1-based 60)."""
    return next(
        seq
        for seq, (parent, name) in RD.mutated_sequences.items()
        if parent == unit and name == "dupC"
    )


def allele(inner: list[str]) -> str:
    """Motif 1..9 sequence: pre-repeats + inner tokens + after-repeats.

    A token is a repeat ID from the bundled dictionary or, if it is not an ID, a raw
    unit sequence (for example ``dupc()``).
    """
    return "".join(RD.repeats.get(u, u) for u in PRE + inner + POST)


def _noisy(seq: str, err: float, rng: random.Random) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < err / 3:
            continue  # deletion
        if r < 2 * err / 3:
            out.append(rng.choice("ACGT"))  # substitution
            continue
        out.append(base)
        if r < err:
            out.append(rng.choice("ACGT"))  # insertion
    return "".join(out)


def reads(
    allele_seq: str,
    n: int,
    *,
    err: float,
    seed: int,
    strand_mix: bool = True,
    flank_bp: int = 40,
    smear_frac: float = 0.0,
    smear_margin_units: int = 10,
    smear_min_deletion_units: int = 5,
) -> list[ReadRecord]:
    """``n`` noisy reads of flank + allele + flank; smear reads lose an internal block.

    The smear cut keeps at least ``smear_margin_units`` repeat units untouched at each
    end and deletes at least ``smear_min_deletion_units`` units; both are capped to a
    fraction of the allele so a short synthetic allele can never make the deletion
    window invalid (the historical fixed-bp margins could, for short alleles).
    """
    rng = random.Random(seed)
    left = RD.flanking_left[-flank_bp:] if flank_bp else ""
    right = RD.flanking_right[:flank_bp] if flank_bp else ""
    unit_bp = RD.repeat_length_bp
    margin = min(smear_margin_units * unit_bp, max(1, len(allele_seq) // 3))
    min_deletion = min(smear_min_deletion_units * unit_bp, max(1, len(allele_seq) // 6))
    out = []
    for i in range(n):
        template = left + allele_seq + right
        if smear_frac and rng.random() < smear_frac:
            cut_a = rng.randint(len(left) + margin, len(left) + len(allele_seq) // 2)
            cut_b = rng.randint(cut_a + min_deletion, len(left) + len(allele_seq) - margin)
            template = template[:cut_a] + template[cut_b:]
        seq = _noisy(template, err, rng)
        if strand_mix and rng.random() < 0.5:
            seq = rc(seq)
        out.append(ReadRecord(f"r{seed}_{i}", seq, "5" * len(seq)))
    return out
