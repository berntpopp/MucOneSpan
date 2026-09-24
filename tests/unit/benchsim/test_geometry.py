"""Case geometry: primers, amplicon intervals, VNTR bounds in both frames."""

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from muc_one_span.benchsim.geometry import (
    amplicon_interval,
    case_geometry,
    flank_sequences,
    haplotype_sequences,
    primer_pair,
    source_inputs,
    vntr_bounds,
)
from muc_one_span.nomenclature import revcomp

FWD, REV = "GGAGAAAAGG", "GCCGTTGTGC"
VNTR1, VNTR2 = "ACACACAC", "ACACACACACAC"
HAPS = {
    1: "tt" + FWD + "CA" + VNTR1 + "AT" + revcomp(REV) + "aa",
    2: "c" + FWD + "CA" + VNTR2 + "AT" + revcomp(REV),
}


def _truth() -> Any:
    return SimpleNamespace(
        haplotypes=(SimpleNamespace(sequence=VNTR1), SimpleNamespace(sequence=VNTR2))
    )


def test_primer_pair_from_config(tmp_path: Path) -> None:
    config = tmp_path / "c.json"
    config.write_text(
        json.dumps({"amplicon_params": {"forward_primer": "acg", "reverse_primer": "tt"}})
    )
    assert primer_pair(config) == ("ACG", "TT")
    config.write_text(json.dumps({"amplicon_params": {"forward_primer": "ACG"}}))
    with pytest.raises(KeyError, match="primers are missing"):
        primer_pair(config)


def test_amplicon_interval_unique_sites() -> None:
    assert amplicon_interval(HAPS[1], FWD, REV) == (2, 2 + 10 + 12 + 10)
    with pytest.raises(ValueError, match="primer sites not unique"):
        amplicon_interval(FWD + FWD + revcomp(REV), FWD, REV)


def test_haplotype_and_flank_sequences(tmp_path: Path) -> None:
    fa = tmp_path / "t.simulated.fa"
    fa.write_text(f">haplotype_2 x\n{HAPS[2]}\n>haplotype_1\n{HAPS[1]}\n")
    assert haplotype_sequences(fa) == {1: HAPS[1].upper(), 2: HAPS[2].upper()}
    flank = tmp_path / "flank.fa"
    flank.write_text(">left\naaa\n>right\nCC\n")
    assert flank_sequences(flank) == ("AAA", "CC") and flank_sequences(None) == ("", "")


def _haps() -> dict[int, str]:
    return {h: s.upper() for h, s in HAPS.items()}


def test_vntr_bounds_validates_truth() -> None:
    haps = {1: "TT" + VNTR1 + "GG", 2: "TT" + VNTR2 + "GG"}
    assert vntr_bounds(_truth(), "TT", haps) == {1: (2, 10), 2: (2, 14)}
    with pytest.raises(ValueError, match="does not match"):
        vntr_bounds(_truth(), "T", haps)


def test_amplicon_geometry_round_trip() -> None:
    haps = _haps()
    vntr = {1: (14, 22), 2: (13, 25)}
    geo = case_geometry("ont_amplicon_r10", haps, vntr, primers=(FWD, REV), flank_ext=(0, 0))
    geo = json.loads(json.dumps(geo))  # as stored in case.json
    assert geo["frame"] == "amplicon" and geo["amplicon"] == {"1": [2, 34], "2": [1, 37]}
    assert geo["vntr_source"] == {"1": [12, 20], "2": [12, 24]} and geo["flank_ext"] is None
    sources, bounds = source_inputs(geo, haps)
    assert sources[1] == FWD + "CA" + VNTR1 + "AT" + revcomp(REV)
    assert sources[2][bounds[2][0] : bounds[2][1]] == VNTR2
    with pytest.raises(ValueError, match="primers are required"):
        case_geometry("hifi_amplicon", haps, vntr, primers=None, flank_ext=(0, 0))


def test_genomic_geometry_needs_matching_flanks() -> None:
    haps = _haps()
    vntr = {1: (14, 22), 2: (13, 25)}
    geo = case_geometry("ont_genomic_targeted", haps, vntr, primers=None, flank_ext=(3, 2))
    assert geo["frame"] == "flanked_source" and geo["primers"] is None
    assert geo["vntr_source"] == {"1": [17, 25], "2": [16, 28]} and geo["flank_ext"] == [3, 2]
    sources, bounds = source_inputs(geo, haps, ("AAA", "CC"))
    assert sources[1] == "AAA" + haps[1] + "CC"
    assert sources[1][bounds[1][0] : bounds[1][1]] == VNTR1
    with pytest.raises(ValueError, match="--flank-fasta"):
        source_inputs(geo, haps)
