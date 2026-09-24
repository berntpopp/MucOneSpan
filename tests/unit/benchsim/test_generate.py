"""Case generation: validation, read truth, manifest and ledger (MucOneUp mocked)."""

import gzip
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, AmountConfig, BenchConfig
from muc_one_span.benchsim.design import Design, build_split
from muc_one_span.benchsim.generate import GenerateContext, generate_case, write_manifest
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation.models import Event, TruthHaplotype, TruthSample

MOD = "muc_one_span.benchsim.generate"
HEADER = "read_id\thap\tmolecule\tkind\tstrand\tsrc_start\tsrc_end\tn_hp_edits\thp_edits\tdetail\n"
FASTQ = {
    "ont_amplicon_r10": "{}_amplicon_ont.fastq",
    "hifi_amplicon": "{}_amplicon_hifi.fastq",
    "ont_genomic_targeted": "{}_ont_fragments.fastq",
}


def _ctx(tmp_path: Path, **kw: Any) -> GenerateContext:
    prof = tmp_path / "profiles"
    prof.mkdir(exist_ok=True)
    for name in ("ont_r10_sup_amplicon_v1", "ont_r10_genomic_v1", "hifi_amplicon_v1"):
        (prof / f"{name}.json").write_text(
            json.dumps({"schema_version": 1, "name": name, "platform": "ont", "molecules": {}})
        )
    config = tmp_path / "c.json"
    primers = {"forward_primer": "TT", "reverse_primer": "CC"}  # revcomp CC = GG flank
    config.write_text(json.dumps({"amplicon_params": primers}))
    ctx = GenerateContext("muconeup", config, tmp_path / "out", prof, None, None, "0.45.0")
    return replace(ctx, **kw)


def _plain(designs: list[Design], profile: str = "ont_amplicon_r10", event: bool = True) -> Design:
    """A Markov design that needs no structure pre-run (see Ruling on 0_identical)."""
    return next(
        d
        for d in designs
        if bool(d.event) == event
        and d.profile == profile
        and d.composition == "markov"
        and d.delta_class != "0_identical"
    )


DEV = build_split("dev", 30, "s", ["dupC"])


def _rd():
    """Tiny dictionary adapted from tests/unit/test_evaluation_truth.py::fixture."""
    return replace(
        load_repeat_dictionary(),
        repeats={"X": "AC"},
        flanking_left="TT",
        flanking_right="GG",
        mutations={
            "dupC": {
                "allowed_repeats": ["X"],
                "changes": [{"type": "insert", "start": 2, "sequence": "C"}],
            }
        },
    )


def _opt(args: list[str], flag: str) -> str:
    return args[args.index(flag) + 1]


def _write_truth(out: Path, base: str, mutant_hap: int | None) -> None:
    seqs = ["ACC" if mutant_hap == i else "AC" for i in (1, 2)]
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{base}.001.simulated.fa").write_text(
        "".join(f">haplotype_{i}\nTT{s}GG\n" for i, s in enumerate(seqs, 1))
    )
    chains = ["Xm" if mutant_hap == i else "X" for i in (1, 2)]
    (out / f"{base}.001.vntr_structure.txt").write_text(
        "".join(f"haplotype_{i}\t{c}\n" for i, c in enumerate(chains, 1))
    )
    rows = [
        {
            "repeat_count": 1,
            "vntr_length": len(s),
            "repeat_lengths": [2],
            "mutant_repeat_count": int(mutant_hap == i),
            "mutation_details": [{"position": 1, "repeat": "X"}] if mutant_hap == i else [],
        }
        for i, s in enumerate(seqs, 1)
    ]
    info = {"mutation_name": "dupC", "mutation_targets": [f"{mutant_hap},1"]} if mutant_hap else {}
    stats = {"haplotype_statistics": rows, "mutation_info": info, "provenance": {"seed": 1}}
    (out / f"{base}.001.simulation_stats.json").write_text(json.dumps(stats))
    if mutant_hap:
        (out / f"{base}.001.mutated_unit.fa").write_text(f">haplotype_{mutant_hap}_repeat_1\nACC\n")


class FakeMucOneUp:
    """Writes minimal MucOneUp outputs; records every call."""

    def __init__(self, fastq_records: int = 4, kind: str = "full") -> None:
        self.calls: list[list[str]] = []
        self.fastq_records = fastq_records
        self.kind = kind

    def __call__(self, args: list[str], **_: Any) -> str:
        self.calls.append(list(args))
        out, base = Path(_opt(args, "--out-dir")), _opt(args, "--out-base")
        if "simulate" in args:
            hap = (
                int(_opt(args, "--mutation-targets").split(",")[0])
                if "--mutation-name" in args
                else None
            )
            _write_truth(out, base, hap)
            return ""
        profile = next(p for p, n in FASTQ.items() if ("amplicon" in args) == ("amplicon" in n))
        if "--platform" in args and _opt(args, "--platform") == "pacbio":
            profile = "hifi_amplicon"
        out.mkdir(parents=True, exist_ok=True)
        (out / FASTQ[profile].format(base)).write_text(
            "".join(f"@r{i}\nAC\n+\nII\n" for i in range(self.fastq_records))
        )
        with gzip.open(out / f"{base}_read_truth.tsv.gz", "wt") as fh:
            fh.write(HEADER)
            for i, hap in enumerate((1, 1, 2, 2), 1):
                fh.write(f"r{i}\t{hap}\t{i}\t{self.kind}\t+\t0\t7\t0\t\t\n")
        (out / f"{base}_amplicon_ont_metadata.tsv").write_text(
            "Parameter\tValue\nCommand\tmuconeup reads --coverage 9 --seed 3\n"
        )
        return ""


def test_invalid_design_is_recorded_not_scored(tmp_path: Path) -> None:
    design = _plain(DEV)
    with patch(f"{MOD}.run_tool"), patch(f"{MOD}.load_truth", side_effect=ValueError("no event")):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "design_invalid" and "no event" in case["error"]
    saved = json.loads((tmp_path / "out/dev" / design.design_id / "case.json").read_text())
    assert saved["status"] == "design_invalid"
    assert not (tmp_path / "out/dev/ledger_sealed.jsonl").exists()


def test_event_mismatch_is_design_invalid(tmp_path: Path) -> None:
    design = _plain(DEV)
    truth = TruthSample(
        "x",
        (
            TruthHaplotype("haplotype_1", "A", ("X",)),
            TruthHaplotype("haplotype_2", "A", ("X",)),
        ),
    )
    with patch(f"{MOD}.run_tool"), patch(f"{MOD}.load_truth", return_value=truth):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "design_invalid" and "event" in case["error"]


def test_event_on_normal_design_is_invalid(tmp_path: Path) -> None:
    design = _plain(DEV, event=False)
    ev = Event(3, "X", "dupC")
    truth = TruthSample("x", (TruthHaplotype("haplotype_1", "A", ("X",), (ev,)),))
    with patch(f"{MOD}.run_tool"), patch(f"{MOD}.load_truth", return_value=truth):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "design_invalid" and "normal" in case["error"]


def test_ok_case_records_depth_hashes_ledger_and_is_idempotent(tmp_path: Path) -> None:
    design = _plain(DEV)
    fake = FakeMucOneUp()
    ctx = _ctx(tmp_path)
    with (
        patch(f"{MOD}.run_tool", side_effect=fake),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        case = generate_case(design, ctx)
        assert case["status"] == "ok", case.get("error")
        assert case["realized_depth"] == {"1": 2, "2": 2}
        assert len(case["profile_sha256"]) == 64
        assert case["muconeup_version"] == "0.45.0"
        assert case["actual_targets"] == [[design.targets[0][0], 1]]
        assert case["composition_effective"] == "markov"
        assert set(case["hashes"]) == {"truth_fa", "fastq", "read_truth"}
        geo = case["geometry"]
        assert geo["frame"] == "amplicon" and geo["primers"] == {"forward": "TT", "reverse": "CC"}
        long_hap = str(design.targets[0][0])  # dupC adds one base on the mutant haplotype
        assert geo["amplicon"][long_hap] == [0, 7] and geo["vntr_haplotype"][long_hap] == [2, 5]
        assert geo["vntr_source"][long_hap] == [2, 5] and geo["flank_ext"] is None
        reads = next(c for c in fake.calls if "reads" in c)
        profile_arg = Path(_opt(reads, "--read-profile"))
        assert profile_arg.parent == ctx.out_root / "profiles"
        split_dir = ctx.out_root / "dev"
        truth_dir = split_dir / design.design_id / "truth"
        assert len(list(truth_dir.glob("*_read_truth.tsv.gz"))) == 1
        assert len(list(truth_dir.glob("*_metadata.tsv"))) == 1
        ledger = (split_dir / "ledger_sealed.jsonl").read_text().splitlines()
        assert len(ledger) == 1
        entry = json.loads(ledger[0])
        assert entry["category"] == f"ont_amplicon_r10:{design.delta_class}:dupC"
        assert entry["token"] == design.design_id
        n_calls = len(fake.calls)
        again = generate_case(design, ctx)
    assert again == case and len(fake.calls) == n_calls


def test_fastq_manifest_mismatch_fails_generation(tmp_path: Path) -> None:
    design = _plain(DEV, event=False)
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp(fastq_records=3)),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "generation_failed" and "records" in case["error"]


def test_amplicon_without_primers_is_generation_failed(tmp_path: Path) -> None:
    design = _plain(DEV, event=False)
    ctx = _ctx(tmp_path)
    ctx.config.write_text("{}")
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp()),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        case = generate_case(design, ctx)
    assert case["status"] == "generation_failed" and "primers" in case["error"]


def test_tool_failure_is_generation_failed(tmp_path: Path) -> None:
    design = _plain(DEV, event=False)
    with patch(f"{MOD}.run_tool", side_effect=RuntimeError("pbsim crashed")):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "generation_failed" and "pbsim crashed" in case["error"]


def test_genomic_case_uses_span_and_n_reads(tmp_path: Path) -> None:
    design = _plain(DEV, "ont_genomic_targeted", event=False)
    flank = tmp_path / "flank.fa"
    flank.write_text(">left\n" + "A" * 50 + "\n>right\n" + "C" * 50 + "\n")
    fake = FakeMucOneUp(kind="fragment")
    with (
        patch(f"{MOD}.run_tool", side_effect=fake),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        case = generate_case(design, _ctx(tmp_path, flank_fasta=flank))
    assert case["status"] == "ok", case.get("error")
    reads = next(c for c in fake.calls if "reads" in c)
    assert int(_opt(reads, "--n-reads")) == case["requested_amount"] > 0
    assert _opt(reads, "--flank-fasta") == str(flank)
    # fragment rows cover [0, 7) of the flanked source; the VNTR sits at [52, 54)
    assert case["span"] == {"1": [52, 54], "2": [52, 54]}
    geo = case["geometry"]
    assert geo["frame"] == "flanked_source" and geo["primers"] is None
    assert geo["vntr_haplotype"]["1"] == [2, 4] and geo["vntr_source"]["1"] == [52, 54]
    assert geo["flank_ext"] == [50, 50] and geo["amplicon"] is None
    assert case["realized_depth"] == {"1": 0, "2": 0}


def test_amplicon_amount_uses_artefacts_and_pcr_share(tmp_path: Path) -> None:
    design = replace(_plain(DEV, event=False), smear=0.25, chimera=0.05, depth=10, pcr="none")
    ctx = _ctx(tmp_path)
    base = ctx.profile_dir / "ont_r10_sup_amplicon_v1.json"
    data = json.loads(base.read_text())
    data["molecules"] = {"concatemer_rate": 0.2}
    base.write_text(json.dumps(data))
    fake = FakeMucOneUp()
    with (
        patch(f"{MOD}.run_tool", side_effect=fake),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        case = generate_case(design, ctx)
    # 10 / (0.5 * (1 - 0.5)) = 40 templates
    assert case["requested_amount"] == 40 and case["amount_capped"] is False
    reads = next(c for c in fake.calls if "reads" in c)
    assert _opt(reads, "--coverage") == "40"


def test_amplicon_amount_is_capped_at_the_minor_share_floor(tmp_path: Path) -> None:
    # Strong PCR bias and a large length difference give a minor share ~1e-4; the
    # uncapped template count (~10^5-10^6) never finishes. The floor bounds it and the
    # minor allele is left below its target depth (recorded, like allelic dropout).
    design = replace(_plain(DEV, event=False), smear=0.25, chimera=0.25, depth=10, pcr="strong")
    fake = FakeMucOneUp()
    with (
        patch(f"{MOD}.run_tool", side_effect=fake),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
        patch(f"{MOD}.pcr_minor_share", return_value=1e-4),
    ):
        case = generate_case(design, _ctx(tmp_path))
    floor = DEFAULT_BENCH_CONFIG.amount.min_minor_share
    assert case["requested_amount"] == math.ceil(10 / (floor * (1 - 0.5)))
    assert case["amount_capped"] is True and case["min_minor_share"] == floor


def test_minor_share_floor_is_configured_and_recorded(tmp_path: Path) -> None:
    design = replace(_plain(DEV, event=False), smear=0.25, chimera=0.25, depth=10, pcr="strong")
    floor = DEFAULT_BENCH_CONFIG.amount.min_minor_share * 2
    cfg = BenchConfig(amount=AmountConfig(min_minor_share=floor))
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp()),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
        patch(f"{MOD}.pcr_minor_share", return_value=1e-4),
    ):
        case = generate_case(design, _ctx(tmp_path, bench=cfg))
    assert case["requested_amount"] == math.ceil(10 / (floor * (1 - 0.5)))
    assert case["min_minor_share"] == floor and case["bench_config_sha256"] == cfg.sha256()


def test_case_records_target_clamp(tmp_path: Path) -> None:
    design = replace(_plain(DEV), target_clamped=True)
    with patch(f"{MOD}.run_tool"), patch(f"{MOD}.load_truth", side_effect=ValueError("x")):
        case = generate_case(design, _ctx(tmp_path))
    assert case["target_clamped"] is True and case["design"]["target_clamped"] is True


def test_pool_structure_too_short_for_an_event_is_design_invalid(tmp_path: Path) -> None:
    pool = tmp_path / "pool.txt"
    pool.write_text("haplotype_1\t1-X-9\nhaplotype_2\t1-X-9\n")
    design = replace(_plain(DEV), composition="real_derived")
    with patch(f"{MOD}.run_tool"):
        case = generate_case(design, _ctx(tmp_path, structure_pool=pool))
    assert case["status"] == "design_invalid" and "too short" in case["error"]


def test_rerun_after_failure_clears_stale_outputs(tmp_path: Path) -> None:
    design = _plain(DEV, event=False)
    ctx = _ctx(tmp_path)
    stale = ctx.out_root / "dev" / design.design_id / "truth"
    stale.mkdir(parents=True)
    (stale / "old.simulated.fa").write_text(">haplotype_1\nA\n")
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp()),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        case = generate_case(design, ctx)
    assert case["status"] == "ok", case.get("error")


def test_write_manifest_sorted(tmp_path: Path) -> None:
    path = write_manifest(tmp_path, [{"design_id": "b"}, {"design_id": "a", "status": "ok"}])
    lines = [json.loads(x) for x in path.read_text().splitlines()]
    assert path.name == "manifest.jsonl" and [x["design_id"] for x in lines] == ["a", "b"]


@pytest.mark.parametrize("composition", ["real_derived"])
def test_real_derived_without_pool_falls_back_to_markov(tmp_path: Path, composition: str) -> None:
    design = replace(_plain(DEV, event=False), composition=composition)
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp()),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "ok", case.get("error")
    assert case["composition_effective"] == "markov"
    assert case["design"]["composition"] == "real_derived"


def test_identical_design_masks_event_positions_only() -> None:
    from muc_one_span.benchsim.generate import DesignInvalidError, validate_events

    design = replace(_plain(DEV), delta_class="0_identical", targets=((2, 2),))
    ev = Event(2, "X", "dupC")
    same = TruthSample(
        "x",
        (
            TruthHaplotype("haplotype_1", "A", ("1", "C", "9")),
            TruthHaplotype("haplotype_2", "A", ("1", "X:dupC", "9"), (ev,)),
        ),
    )
    assert validate_events(same, design) == [[2, 2]]
    differ = replace(
        same,
        haplotypes=(replace(same.haplotypes[0], structure=("1", "C", "8")), same.haplotypes[1]),
    )
    with pytest.raises(DesignInvalidError, match="differ"):
        validate_events(differ, design)


def test_identical_normal_design_checks_structures() -> None:
    from muc_one_span.benchsim.generate import DesignInvalidError, validate_events

    design = replace(_plain(DEV, event=False), delta_class="0_identical")
    haps = (
        TruthHaplotype("haplotype_1", "A", ("1", "X", "9")),
        TruthHaplotype("haplotype_2", "A", ("1", "A", "9")),
    )
    with pytest.raises(DesignInvalidError, match="differ"):
        validate_events(TruthSample("x", haps), design)
    assert validate_events(TruthSample("x", (haps[0], haps[0])), design) == []


def test_resume_refuses_a_case_made_under_other_settings(tmp_path: Path) -> None:
    from muc_one_span.benchsim.generate import StaleCaseError

    design = _plain(DEV, event=False)
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp()),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        first = generate_case(design, _ctx(tmp_path))
        assert first["status"] == "ok", first.get("error")
        assert generate_case(design, _ctx(tmp_path)) == first  # same settings: reused
        other = BenchConfig(amount=AmountConfig(genomic_mc_draws=7))
        with pytest.raises(StaleCaseError, match="--out-root"):
            generate_case(design, _ctx(tmp_path, bench=other))
        with pytest.raises(StaleCaseError, match="design"):
            generate_case(replace(design, depth=design.depth + 1), _ctx(tmp_path))


def test_resume_hashes_only_generation_settings(tmp_path: Path) -> None:
    from muc_one_span.benchsim.bench_config import AtlasConfig, DesignConfig, ReportConfig
    from muc_one_span.benchsim.generate import StaleCaseError

    design = _plain(DEV, event=False)
    base = DEFAULT_BENCH_CONFIG
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp()),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        first = generate_case(design, _ctx(tmp_path))
        assert first["bench_generation_sha256"] == base.generation_sha256()
        assert first["bench_config_sha256"] == base.sha256()
        atlas_only = replace(
            base,
            atlas=AtlasConfig(min_resolvable_depth=base.atlas.min_resolvable_depth + 1),
            report=ReportConfig(bootstrap_seed=base.report.bootstrap_seed + 1),
        )
        assert generate_case(design, _ctx(tmp_path, bench=atlas_only)) == first
        levels = base.design.chimera_levels[:1]
        design_cfg = replace(base, design=DesignConfig(chimera_levels=levels))
        with pytest.raises(StaleCaseError, match="generation settings"):
            generate_case(design, _ctx(tmp_path, bench=design_cfg))


def test_resume_accepts_a_legacy_case_with_the_same_full_hash(tmp_path: Path) -> None:
    from muc_one_span.benchsim.generate import StaleCaseError

    design = _plain(DEV, event=False)
    with (
        patch(f"{MOD}.run_tool", side_effect=FakeMucOneUp()),
        patch(f"{MOD}.load_repeat_dictionary", return_value=_rd()),
    ):
        first = generate_case(design, _ctx(tmp_path))
        case_json = tmp_path / "out" / design.split / design.design_id / "case.json"
        legacy = {k: v for k, v in first.items() if k != "bench_generation_sha256"}
        case_json.write_text(json.dumps(legacy))
        assert generate_case(design, _ctx(tmp_path)) == legacy
        other = BenchConfig(atlas=replace(DEFAULT_BENCH_CONFIG.atlas, top_reasons=1))
        with pytest.raises(StaleCaseError, match="generation settings"):
            generate_case(design, _ctx(tmp_path, bench=other))
