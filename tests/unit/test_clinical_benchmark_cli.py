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
