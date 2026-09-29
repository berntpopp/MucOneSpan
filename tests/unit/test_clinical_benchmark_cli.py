"""Portable clinical benchmark command contracts."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def test_cli_requires_subcommand():
    import pytest

    spec = spec_from_file_location("clinical_benchmark", Path("scripts/clinical_benchmark.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    with pytest.raises(SystemExit) as error:
        module.main([])
    assert error.value.code == 2


def test_preparation_error_does_not_skip_later_runs(tmp_path, monkeypatch):
    import json

    from muc_one_span.clinical_data import inventory_from_tsv

    spec = spec_from_file_location("clinical_benchmark", Path("scripts/clinical_benchmark.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    header = "run_accession\tsample_accession\texperiment_accession\tsample_alias\tinstrument_platform\tlibrary_strategy\tfastq_ftp\tfastq_bytes\tfastq_md5\n"
    rows = "".join(
        f"ERR{i}\tSAMEA{i}\tERX{i}\tHG00{i}_PCR_MUC1\tOXFORD_NANOPORE\tAMPLICON\tftp.example/{i}.gz\t1\t"
        + "0" * 32
        + "\n"
        for i in (1, 2)
    )
    inventory = inventory_from_tsv(
        header + rows, study="PRJEB1", source_url="https://example", retrieved_at="2026-09-15"
    )
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(inventory))
    attempted = []

    def prepare(run, root):
        attempted.append(run["run_accession"])
        if run["run_accession"] == "ERR1":
            raise ValueError("synthetic corrupt input")
        (root / run["run_accession"]).mkdir(parents=True)
        return {"run_accession": run["run_accession"]}

    monkeypatch.setattr(module, "prepare_run", prepare)
    assert (
        module.main(["prepare", "--manifest", str(manifest), "--data-root", str(tmp_path / "data")])
        == 1
    )
    assert attempted == ["ERR1", "ERR2"]
    assert (tmp_path / "data/ERR2/preparation.json").is_file()
    records = json.loads((tmp_path / "data/cohort_attempts.json").read_text())["records"]
    assert records[0]["status"] == "execution_failed"
    assert records[0]["error"] == "synthetic corrupt input"


def test_run_defaults_to_hybrid_and_keeps_legacy_ladder_hash(tmp_path, monkeypatch):
    import json

    spec = spec_from_file_location("clinical_benchmark", Path("scripts/clinical_benchmark.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"runs": []}))
    environment = tmp_path / "environment.json"
    environment.write_text(json.dumps({"model": {"path": "model.path"}}))
    data_root = tmp_path / "data"
    (data_root / "ERR1").mkdir(parents=True)
    (data_root / "ERR1" / "preparation.json").write_text(json.dumps({"run_accession": "ERR1"}))

    run = {"run_accession": "ERR1", "arm": "primary_amplicon"}
    monkeypatch.setattr(module, "validate_inventory", lambda raw: {"runs": [run]})
    monkeypatch.setattr(module, "verify_environment", lambda environment, checkout: None)

    captured: list[dict] = []

    def fake_run_case(run_arg, preparation, output_root, settings, *, resume=False):
        captured.append(settings)
        return {
            "run_accession": run_arg["run_accession"],
            "status": "completed",
            "analysis_state": "completed",
        }

    monkeypatch.setattr(module, "run_case", fake_run_case)

    def invoke(output_root, extra):
        args = [
            "run",
            "--manifest",
            str(manifest),
            "--data-root",
            str(data_root),
            "--output-root",
            str(output_root),
            "--environment",
            str(environment),
            *extra,
        ]
        assert module.main(args) == 0

    invoke(tmp_path / "out_default", [])
    assert captured[-1]["engine"] == "hybrid"
    assert "assay" not in captured[-1]

    # The ladder hashes without an engine key, so pre-0.17 ladder attempts stay resumable.
    invoke(tmp_path / "out_ladder", ["--engine", "ladder"])
    assert "engine" not in captured[-1]

    invoke(tmp_path / "out_hybrid", ["--engine", "hybrid", "--assay", "genomic"])
    assert captured[-1]["engine"] == "hybrid"
    assert captured[-1]["assay"] == "genomic"


def test_load_records_merges_terminal_files_and_latest_failure_journal(tmp_path):
    import json

    spec = spec_from_file_location("clinical_benchmark", Path("scripts/clinical_benchmark.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    for accession in ["ERR1", "ERR2"]:
        root = tmp_path / accession
        root.mkdir()
        (root / "benchmark_record.json").write_text(
            json.dumps({"run_accession": accession, "status": "completed"})
        )
    (tmp_path / "cohort_attempts.json").write_text(
        json.dumps({"records": [{"run_accession": "ERR1", "status": "execution_failed"}]})
    )
    records = {r["run_accession"]: r for r in module.load_records(tmp_path)}
    assert records["ERR1"]["status"] == "execution_failed"
    assert records["ERR2"]["status"] == "completed"


def test_freeze_model_is_optional_and_a_modelless_environment_runs_hybrid_only(
    tmp_path, monkeypatch
):
    import json

    import pytest

    from muc_one_span import clinical_provenance

    spec = spec_from_file_location("clinical_benchmark", Path("scripts/clinical_benchmark.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    args = module.parser().parse_args(["freeze", "--output", str(tmp_path / "e.json")])
    assert args.model is None

    frozen: list[object] = []
    monkeypatch.setattr(
        module, "freeze_environment", lambda checkout, model: frozen.append(model) or {}
    )
    assert module.main(["freeze", "--output", str(tmp_path / "e.json")]) == 0
    assert frozen == [None]

    # verify_environment accepts a model-less environment (no model to compare).
    env = {"model": None}
    assert clinical_provenance.frozen_model_path(env) is None

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"runs": []}))
    environment = tmp_path / "environment.json"
    environment.write_text(json.dumps(env))
    data_root = tmp_path / "data"
    (data_root / "ERR1").mkdir(parents=True)
    (data_root / "ERR1" / "preparation.json").write_text(json.dumps({"run_accession": "ERR1"}))
    run = {"run_accession": "ERR1", "arm": "primary_amplicon"}
    monkeypatch.setattr(module, "validate_inventory", lambda raw: {"runs": [run]})
    monkeypatch.setattr(module, "verify_environment", lambda environment, checkout: None)
    captured: list[dict] = []

    def fake_run_case(run_arg, preparation, output_root, settings, *, resume=False):
        captured.append(settings)
        return {"run_accession": "ERR1", "status": "completed", "analysis_state": "completed"}

    monkeypatch.setattr(module, "run_case", fake_run_case)
    base = ["run", "--manifest", str(manifest), "--data-root", str(data_root)]
    base += ["--environment", str(environment)]
    assert module.main([*base, "--output-root", str(tmp_path / "h")]) == 0
    assert captured[-1]["model"] is None
    with pytest.raises(SystemExit, match="model"):
        module.main([*base, "--output-root", str(tmp_path / "l"), "--engine", "ladder"])
