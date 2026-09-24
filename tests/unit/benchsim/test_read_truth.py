"""Tests for the MucOneUp per-read truth manifest loader."""

import gzip
from pathlib import Path

import pytest

from muc_one_span.benchsim.read_truth import composition, load_read_truth, realized_depth

HEADER = "read_id\thap\tmolecule\tkind\tstrand\tsrc_start\tsrc_end\tn_hp_edits\thp_edits\tdetail\n"


def _write(path: Path, rows: list[str], header: str = HEADER) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write(header + "".join(r + "\n" for r in rows))
    return path


def test_load_and_amplicon_depth(tmp_path: Path) -> None:
    p = _write(
        tmp_path / "t.tsv.gz",
        [
            "b_h1_m0000001\t1\t1\tfull\t+\t0\t2500\t0\t\t",
            "b_h2_m0000002\t2\t2\tsmear\t-\t0\t2500\t0\t\tdeletion:100-900",
            "b_h2_m0000003\t2\t3\tfull\t-\t0\t2600\t1\t52:C:7>8\t",
            "b_h1_m0000004\t1\t4\tofftarget\t+\t10\t300\t0\t\t",
        ],
    )
    rows = load_read_truth(p)
    assert realized_depth(rows, None) == {1: 1, 2: 1}
    comp = composition(rows)
    assert comp["kind_frac"]["offtarget"] == 0.25
    assert comp["strand_frac"]["-"] == 0.5


def test_fragment_depth_requires_covering_span(tmp_path: Path) -> None:
    p = _write(
        tmp_path / "t.tsv.gz",
        [
            "b_h1_m0000001\t1\t1\tfragment\t+\t900\t7000\t0\t\t",
            "b_h1_m0000002\t1\t2\tfragment\t+\t1100\t7000\t0\t\t",
        ],
    )
    assert realized_depth(load_read_truth(p), {1: (1000, 5000), 2: (1000, 6000)}) == {1: 1, 2: 0}


def test_bad_header_and_duplicates(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="header"):
        load_read_truth(_write(tmp_path / "a.tsv.gz", [], header="read_id\thap\n"))
    row = "b_h1_m0000001\t1\t1\tfull\t+\t0\t10\t0\t\t"
    with pytest.raises(ValueError, match="duplicate"):
        load_read_truth(_write(tmp_path / "b.tsv.gz", [row, row]))
