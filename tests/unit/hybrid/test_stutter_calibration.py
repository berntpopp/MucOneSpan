"""Task 15f fix round 1: pinned pre-15f numbers and the tolerated wild-type share.

``test_shift_model_reproduces_pre_15f_numbers``: the expected read support was computed
with the evidence code of 9f39eff (the commit before Task 15f) on these exact fixtures
(``test_evidence_stutter._reads``, seed ``SEEDS[0]``) and committed as literals, so
``hp_stutter_model = "shift"`` provably is the former model.

``test_tolerated_wild_type_share``: at ``event_max_alternative_frac`` 0.25, how many of
the three seeds still rate a wild-type/dupC mixture ``supported``, per stutter shape and
model. It documents the tolerated wild-type share for calibration; a change to the
mixture fit that moves any cell must update this table on purpose.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable
from typing import Any

import pytest

from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import test_evidence_stutter as base
from tests.unit.hybrid import test_stutter_guard as guard

S = HybridSettings()
SHIFT = dataclasses.replace(S, hp_stutter_model="shift")
FIELDS = ("n", "alt", "ref", "other", "alt_frac", "alternative_frac", "llr", "strand_llr")
PRE_15F: dict[str, tuple[list[tuple[str, int]], Callable[[int], float], dict[str, Any]]] = {
    "log_linear_dupc": (
        [(base.DUPC, base.N_READS)],
        base._log_linear,
        {"n": 200, "alt": 117, "ref": 65, "other": 18, "alt_frac": 0.585,
         "alternative_frac": 0.23, "llr": 234.4, "strand_llr": {"+": 123.1, "-": 111.2},
         "status": "supported"},
    ),
    "d1_step_dupc": (
        [(base.DUPC, base.N_READS)],
        base._d1_step,
        {"n": 200, "alt": 118, "ref": 67, "other": 15, "alt_frac": 0.59,
         "alternative_frac": 0.287, "llr": 197.8, "strand_llr": {"+": 111.4, "-": 86.4},
         "status": "discordant"},
    ),
    "wild_type_reads": (
        [(base.WILD, base.N_READS)],
        base._log_linear,
        {"n": 200, "alt": 15, "ref": 166, "other": 19, "alt_frac": 0.075,
         "alternative_frac": 0.961, "llr": -268.9, "strand_llr": {"+": -155.6, "-": -113.3},
         "status": "not_supported"},
    ),
    "wild_type_40pct": (
        [(base.DUPC, 120), (base.WILD, 80)],
        base._log_linear,
        {"n": 200, "alt": 81, "ref": 100, "other": 19, "alt_frac": 0.405,
         "alternative_frac": 0.489, "llr": 58.0, "strand_llr": {"+": 11.9, "-": 46.1},
         "status": "discordant"},
    ),
}  # fmt: skip


@pytest.mark.parametrize("case", sorted(PRE_15F))
def test_shift_model_reproduces_pre_15f_numbers(case: str) -> None:
    parts, p_del, expected = PRE_15F[case]
    reads = base._reads(parts, base.SEEDS[0], p_del)
    support = base._dupc_support(base.DUPC, reads, SHIFT)[1][0]["read_support"]
    assert {k: support[k] for k in (*FIELDS, "status")} == expected


SHAPES: dict[str, Callable[[int], float]] = {
    "log_linear": base._log_linear,
    "d1_step": base._d1_step,
    "ont_plus": guard._ont_plus,
}
WILD_SHARES = (0.0, 0.2, 0.25, 0.3)
# Seeds (of 3) rated supported, per (shape, model): one count per WILD_SHARES entry.
SUPPORTED_SEEDS = {
    ("log_linear", "length"): (3, 1, 0, 0),
    ("log_linear", "shift"): (3, 0, 0, 0),
    ("d1_step", "length"): (3, 0, 0, 0),
    ("d1_step", "shift"): (2, 0, 0, 0),
    # The saturating ONT "+" shape falls back to the shift model (identifiability
    # guard); both then tolerate up to a 25% wild-type share.
    ("ont_plus", "length"): (3, 3, 3, 0),
    ("ont_plus", "shift"): (3, 3, 3, 0),
}


@pytest.mark.parametrize(("shape", "model"), sorted(SUPPORTED_SEEDS))
def test_tolerated_wild_type_share(shape: str, model: str) -> None:
    s = dataclasses.replace(S, hp_stutter_model=model)
    counts = []
    for share in WILD_SHARES:
        supported = 0
        for seed in base.SEEDS:
            support = guard._support(guard._mix(share, seed, SHAPES[shape]), s)
            supported += support["status"] == "supported"
        counts.append(supported)
    assert tuple(counts) == SUPPORTED_SEEDS[(shape, model)]
