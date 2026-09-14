"""Configuration and fail-closed generation contracts without external tools."""

import json
from pathlib import Path

import pytest

from muc_one_span import experiments as exp


def case(**overrides):
    return {
        "sample": "sample",
        "platform": "hifi",
        "seed": 17,
        "lengths": [25, 30],
        "mutation": "dupC",
        "targets": [[1, 10]],
        "requested_templates": 12,
    } | overrides


def design(tmp_path, cases=None, **options):
    path = tmp_path / "design.json"
    path.write_text(json.dumps({"schema_version": 1, "cases": cases or [case()], **options}))
    return path


def test_explicit_commands_and_platform_configs(tmp_path):
    path = design(
        tmp_path,
        [case(), case(sample="ont", platform="ont")],
        platform_configs={"hifi": "hifi.json", "ont": "ont.json"},
    )
    plan = exp.load_design(path)
    assert plan["platform_configs"]["ont"] == str(tmp_path / "ont.json")
    cmds = exp.case_commands(plan["cases"][1], tmp_path / "ont.json", tmp_path / "out", "sim")
    assert cmds[0].count("--fixed-lengths") == 2
    assert cmds[0][-2:] == ["--mutation-targets", "1,10"]
    assert cmds[1][-2:] == ["--platform", "ont"]
    assert "--track-read-source" not in cmds[1]


@pytest.mark.parametrize(
    "changes",
    [
        {"platform": "pacbio"},
        {"seed": True},
        {"seed": -1},
        {"lengths": [0, 30]},
        {"lengths": [25]},
        {"lengths": [25, False]},
        {"requested_templates": 0},
        {"requested_templates": 2.5},
        {"sample": "../x"},
        {"sample": "a/b"},
        {"sample": "a b"},
        {"targets": [[1, 26]]},
        {"targets": [[3, 1]]},
        {"targets": [[True, 1]]},
        {"targets": [[1, 10], [1, 10]]},
        {"mutation": None},
        {"mutation": ""},
        {"extra": 1},
    ],
)
def test_reject_case_schema(tmp_path, changes):
    with pytest.raises(ValueError):
        exp.load_design(design(tmp_path, [case(**changes)]))


@pytest.mark.parametrize(
    "text",
    [
        '{"schema_version":1,"schema_version":1,"cases":[]}',
        '{"schema_version":true,"cases":[]}',
        '{"schema_version":1,"cases":[],"x":NaN}',
        '{"schema_version":1,"cases":[],"x":1}',
        "[]",
    ],
)
def test_reject_json(tmp_path, text):
    path = tmp_path / "bad.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        exp.load_design(path)


def test_duplicate_names(tmp_path):
    with pytest.raises(ValueError):
        exp.load_design(design(tmp_path, [case(), case()]))


def test_dry_run_has_no_outputs_or_execution(tmp_path, monkeypatch):
    monkeypatch.setattr(exp, "run_tool", lambda *a, **k: pytest.fail("external execution"))
    out = tmp_path / "out"
    result = exp.run_experiment(
        design(tmp_path), out, config=tmp_path / "config.json", dry_run=True
    )
    assert result["cases"][0]["status"] == "planned"
    assert not out.exists()


def test_legacy_panel_preserved():
    plan = exp.load_design(Path("examples/development-experiment.json"))
    assert len(plan["cases"]) == 29
    assert sum(c["platform"] == "hifi" for c in plan["cases"]) == 26
    assert (
        next(c for c in plan["cases"] if c["sample"].endswith("cov50"))["requested_templates"] == 50
    )
    assert [c["seed"] for c in plan["cases"][-3:]] == [1001, 1006, 1002]


def config_file(tmp_path):
    path = tmp_path / "config.json"
    model = tmp_path / "model"
    model.write_text("model")
    path.write_text(
        json.dumps(
            {
                "pacbio_params": {"model_file": "model"},
                "ont_amplicon_params": {"model_file": "model"},
            }
        )
    )
    return path


def test_all_failures_stay_in_manifest(tmp_path, monkeypatch):
    def tool(cmd, **kwargs):
        if "--version" in cmd:
            return "MucOneUp, version test"
        raise RuntimeError("exit code 1: failed")

    monkeypatch.setattr(exp, "run_tool", tool)
    out = tmp_path / "out"
    result = exp.run_experiment(
        design(tmp_path, [case(), case(sample="other")]), out, config=config_file(tmp_path)
    )
    assert [c["status"] for c in result["cases"]] == ["execution_failed"] * 2
    assert len(json.loads((out / "inventory.json").read_text())) == 2
    assert len(result["cases"][0]["commands"]) == 1
    assert result["cases"][0]["commands"][0]["status"] == "execution_failed"
    with pytest.raises(FileExistsError):
        exp.run_experiment(design(tmp_path), out, config=config_file(tmp_path))


def test_fastq_validation_counts_records_not_names(tmp_path):
    p = tmp_path / "reads.fastq"
    p.write_text("@same\nAC\n+\nII\n@same\nTG\n+\nII\n")
    assert exp.count_usable_records(p) == 2
    p.write_text("@bad\nAC\n+\nI\n")
    with pytest.raises(ValueError):
        exp.count_usable_records(p)


def test_bam_validation_uses_primary_records_with_sequence(tmp_path, monkeypatch):
    p = tmp_path / "reads.bam"
    p.touch()
    monkeypatch.setattr(exp, "run_tool", lambda cmd, **kw: "")
    monkeypatch.setattr(
        exp,
        "run_tool_iter",
        lambda cmd, **kw: iter(
            ["same\t0\tx\t1\t60\t2M\t*\t0\t0\tAC\tII", "same\t4\t*\t0\t0\t*\t*\t0\t0\tTG\tII"]
        ),
    )
    assert exp.count_usable_records(p) == 2


def fake_simulator(cmd, **kwargs):
    """Write independently reconstructed diploid simulator-format truth."""
    from muc_one_span.config import _apply_mutation, load_repeat_dictionary

    if "--version" in cmd:
        return "MucOneUp, version fixture"
    out = Path(cmd[cmd.index("--out-dir") + 1])
    if "simulate" not in cmd:
        (out / "reads.fastq").write_text("@same\nAC\n+\nII\n@same\nTG\n+\nII\n")
        return "reads complete"
    rd = load_repeat_dictionary()
    stats, fasta, structures = [], [], []
    mutation = "--mutation-name" in cmd
    for hap, length in ((1, 25), (2, 30)):
        sequences = [rd.repeats["X"]] * length
        labels = ["X"] * length
        details = []
        if mutation and hap == 1:
            sequences[9] = _apply_mutation(rd.repeats["X"], rd.mutations["dupC"]["changes"])
            labels[9] = "Xm"
            details = [{"position": 10, "repeat": "X"}]
            (out / "s.mutated_unit.fa").write_text(f">haplotype_1_repeat_10\n{sequences[9]}\n")
        fasta.append(
            f">haplotype_{hap}\n{rd.flanking_left}{''.join(sequences)}{rd.flanking_right}\n"
        )
        structures.append(f"haplotype_{hap}\t{'-'.join(labels)}\n")
        stats.append(
            {
                "repeat_count": length,
                "vntr_length": sum(map(len, sequences)),
                "repeat_lengths": list(map(len, sequences)),
                "mutant_repeat_count": len(details),
                "mutation_details": details,
            }
        )
    name = cmd[cmd.index("--out-base") + 1]
    (out / f"{name}.001.simulated.fa").write_text("".join(fasta))
    (out / "s.vntr_structure.txt").write_text("".join(structures))
    (out / "s.simulation_stats.json").write_text(
        json.dumps(
            {
                "haplotype_statistics": stats,
                "mutation_info": {"mutation_name": "dupC", "mutation_targets": ["1,10"]}
                if mutation
                else {},
                "provenance": {"seed": 17},
            }
        )
    )
    return "truth complete"


@pytest.mark.parametrize("platform", ["hifi", "ont"])
@pytest.mark.parametrize("mutant", [True, False])
def test_generated_truth_and_actual_inventory(tmp_path, monkeypatch, platform, mutant):
    monkeypatch.setattr(exp, "run_tool", fake_simulator)
    c = case(platform=platform) if mutant else case(platform=platform, mutation=None, targets=[])
    out = tmp_path / "out"
    result = exp.run_experiment(design(tmp_path, [c]), out, config=config_file(tmp_path))
    row = result["cases"][0]
    assert row["status"] == "completed", row.get("error")
    assert row["usable_records"] == 2
    assert row["requested_templates"] == 12
    assert row["config"]["artifact_hashes"][str(tmp_path / "model")]
    inventory = json.loads((out / "inventory.json").read_text())
    assert inventory[0]["input"] == str(out / "sample" / "reads.fastq")
    assert inventory[0]["truth_dir"] == str(out / "sample")


@pytest.mark.parametrize(
    "options",
    [
        {"platform_configs": []},
        {"platform_configs": {"pacbio": "x"}},
        {"platform_configs": {"hifi": None}},
        {"repeat_dictionary": False},
    ],
)
def test_invalid_design_paths(tmp_path, options):
    with pytest.raises(ValueError):
        exp.load_design(design(tmp_path, **options))


def test_config_precedence_and_missing_config(tmp_path):
    p = design(tmp_path, platform_configs={"hifi": "design-config.json"})
    fallback, explicit = tmp_path / "env.json", tmp_path / "cli.json"

    def selected(**kw):
        return exp.run_experiment(p, tmp_path / "out", dry_run=True, **kw)["cases"][0][
            "planned_commands"
        ][0][2]

    assert selected(fallback_config=fallback) == str(tmp_path / "design-config.json")
    assert selected(config=explicit, fallback_config=fallback) == str(explicit)
    p = design(tmp_path)
    assert selected(fallback_config=fallback) == str(fallback)
    with pytest.raises(ValueError, match="no MucOneUp config"):
        selected()


def test_version_failure_is_recorded_before_case_execution(tmp_path, monkeypatch):
    monkeypatch.setattr(
        exp,
        "run_tool",
        lambda *a, **kw: (_ for _ in ()).throw(FileNotFoundError("missing simulator")),
    )
    out = tmp_path / "out"
    report = exp.run_experiment(design(tmp_path), out, config=tmp_path / "x")
    assert report["cases"][0]["status"] == "execution_failed"
    assert report["version_command"] == [report["cases"][0]["planned_commands"][0][0], "--version"]
    assert report["version_result"]["status"] == "execution_failed"


def test_relative_executable_stays_bound_to_invocation_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    report = exp.run_experiment(
        design(tmp_path),
        tmp_path / "out",
        config=tmp_path / "config",
        executable="bin/sim",
        dry_run=True,
    )
    assert report["cases"][0]["planned_commands"][0][0] == str(tmp_path / "bin/sim")


def test_mutated_generated_truth_failure_cannot_be_completed(tmp_path, monkeypatch):
    def tool(cmd, **kwargs):
        output = fake_simulator(cmd, **kwargs)
        if "reads" in cmd:
            p = Path(cmd[cmd.index("--out-dir") + 1]) / "s.vntr_structure.txt"
            p.write_text(p.read_text().replace("Xm", "X"))
        return output

    monkeypatch.setattr(exp, "run_tool", tool)
    out = tmp_path / "out"
    report = exp.run_experiment(design(tmp_path), out, config=config_file(tmp_path))
    assert report["cases"][0]["status"] == "execution_failed"
    assert "usable_records" not in report["cases"][0]
    assert "__generation_failed__" in json.loads((out / "inventory.json").read_text())[0]["input"]


def test_malformed_simulator_section_is_retained_case_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(exp, "run_tool", fake_simulator)
    cfg = config_file(tmp_path)
    cfg.write_text('{"pacbio_params":[]}')
    report = exp.run_experiment(design(tmp_path), tmp_path / "out", config=cfg)
    assert report["cases"][0]["status"] == "execution_failed"
    assert "pacbio_params" in report["cases"][0]["error"]


def test_model_changes_stop_before_reads_and_stay_visible(tmp_path, monkeypatch):
    def tool(cmd, **kwargs):
        result = fake_simulator(cmd, **kwargs)
        if "simulate" in cmd:
            (tmp_path / "model").write_text("changed model")
        return result

    monkeypatch.setattr(exp, "run_tool", tool)
    report = exp.run_experiment(design(tmp_path), tmp_path / "out", config=config_file(tmp_path))
    row = report["cases"][0]
    assert row["status"] == "execution_failed"
    assert len(row["commands"]) == 1
    assert "artifact changed" in row["error"]


def test_model_changed_by_final_reads_cannot_be_completed(tmp_path, monkeypatch):
    def tool(cmd, **kwargs):
        result = fake_simulator(cmd, **kwargs)
        if "reads" in cmd:
            (tmp_path / "model").write_text("changed during final command")
        return result

    monkeypatch.setattr(exp, "run_tool", tool)
    out = tmp_path / "out"
    report = exp.run_experiment(design(tmp_path), out, config=config_file(tmp_path))
    row = report["cases"][0]
    assert row["status"] == "execution_failed"
    assert "artifact changed" in row["error"]
    assert "usable_records" not in row
    assert "__generation_failed__" in json.loads((out / "inventory.json").read_text())[0]["input"]


def test_clean_external_path_pins_version_commands_and_executable_hash(tmp_path, monkeypatch):
    import hashlib

    virtual_bin, external_bin = tmp_path / ".venv/bin", tmp_path / "external/bin"
    for directory, text in ((virtual_bin, "virtualenv program"), (external_bin, "actual program")):
        directory.mkdir(parents=True)
        program = directory / "muconeup"
        program.write_text(text)
        program.chmod(0o755)
    monkeypatch.setenv("PATH", f"{virtual_bin}:{external_bin}")
    commands = []

    def tool(cmd, **kwargs):
        commands.append(cmd)
        return fake_simulator(cmd, **kwargs)

    monkeypatch.setattr(exp, "run_tool", tool)
    report = exp.run_experiment(design(tmp_path), tmp_path / "out", config=config_file(tmp_path))
    actual = str(external_bin / "muconeup")
    assert report["cases"][0]["status"] == "completed"
    assert all(cmd[0] == actual for cmd in commands)
    assert report["simulator_executable"] == actual
    assert report["simulator_executable_sha256"] == hashlib.sha256(b"actual program").hexdigest()
    assert report["version_command"][0] == actual
    assert all(cmd[0] == actual for cmd in report["cases"][0]["planned_commands"])
