from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.benchsim.muconeup import reads_args, require_muconeup, simulate_args

MOD = "muc_one_span.benchsim.muconeup"


@pytest.mark.parametrize(
    "out,ok",
    [
        ("muconeup, version 0.45.0\n", True),
        ("muconeup, version 0.44.5\n", False),
        ("muconeup, version 1.2.3\n", True),
    ],
)
def test_version_gate(out: str, ok: bool) -> None:
    with patch(f"{MOD}.run_tool", return_value=out):
        if ok:
            assert require_muconeup("muconeup") == out.split()[-1]
        else:
            with pytest.raises(RuntimeError, match=r"0\.45\.0"):
                require_muconeup("muconeup")


def test_simulate_args_fixed_lengths_and_mutation(tmp_path: Path) -> None:
    args = simulate_args(
        "mu", tmp_path / "c.json", tmp_path, "b", 7, (40, 61), None, "dupC", ((2, 30),)
    )
    assert args[:4] == ["mu", "--config", str(tmp_path / "c.json"), "simulate"]
    assert args.count("--fixed-lengths") == 2 and "--output-structure" in args
    assert args[args.index("--mutation-name") + 1] == "dupC"
    assert args[args.index("--mutation-targets") + 1] == "2,30"


def test_simulate_args_structure_file_excludes_lengths(tmp_path: Path) -> None:
    s = tmp_path / "s.txt"
    args = simulate_args("mu", tmp_path / "c.json", tmp_path, "b", 7, None, s, None, ())
    assert "--fixed-lengths" not in args and args[args.index("--input-structure") + 1] == str(s)


def test_simulate_args_requires_lengths_or_structure(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="either lengths or structure_file is required"):
        simulate_args("mu", tmp_path / "c.json", tmp_path, "b", 7, None, None, None, ())


def test_reads_args_per_profile(tmp_path: Path) -> None:
    fa = tmp_path / "t.fa"
    amp = reads_args(
        "mu",
        tmp_path / "c.json",
        "ont_amplicon_r10",
        fa,
        tmp_path,
        "b",
        3,
        amount=900,
        profile_path=None,
        flank_fasta=None,
        pcr_preset=None,
    )
    assert amp[3:5] == ["reads", "amplicon"] and "--no-align" in amp
    assert amp[amp.index("--read-profile") + 1] == "ont_r10_sup_amplicon_v1"
    assert amp[amp.index("--coverage") + 1] == "900"
    gen = reads_args(
        "mu",
        tmp_path / "c.json",
        "ont_genomic_targeted",
        fa,
        tmp_path,
        "b",
        3,
        amount=250,
        profile_path=tmp_path / "p.json",
        flank_fasta=tmp_path / "f.fa",
        pcr_preset=None,
    )
    assert gen[3:5] == ["reads", "ont"]
    assert gen[gen.index("--simulator") + 1] == "pbsim3-fragments"
    assert gen[gen.index("--read-profile") + 1] == str(tmp_path / "p.json")
    assert gen[gen.index("--n-reads") + 1] == "250"
    assert gen[gen.index("--flank-fasta") + 1] == str(tmp_path / "f.fa")
    hifi = reads_args(
        "mu",
        tmp_path / "c.json",
        "hifi_amplicon",
        fa,
        tmp_path,
        "b",
        3,
        amount=60,
        profile_path=None,
        flank_fasta=None,
        pcr_preset="no_bias",
    )
    assert hifi[hifi.index("--read-profile") + 1] == "hifi_amplicon_v1"
    assert hifi[hifi.index("--pcr-preset") + 1] == "no_bias"
