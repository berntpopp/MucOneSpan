"""Clinical runner lifecycle and provenance regression tests."""

from pathlib import Path

import pytest

from muc_one_span.clinical_runner import run_case, validate_record


def test_rejects_unverifiable_completed_record():
    assert not validate_record({"status": "completed", "events": [], "outputs": {}})


def test_rejects_changed_input_before_execution(tmp_path: Path):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="input"):
        run_case(
            {"run_accession": "ERR1", "arm": "primary_amplicon"},
            {"run_accession": "ERR1", "output_path": str(reads), "output_sha256": "stale"},
            tmp_path / "results",
            {"threads": 2, "timeout": 10},
        )


def test_refuses_nonempty_attempt_directory(tmp_path: Path):
    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    result = tmp_path / "results" / "ERR1"
    result.mkdir(parents=True)
    (result / "summary.json").write_text("{}")
    with pytest.raises(ValueError, match="existing"):
        run_case(
            {"run_accession": "ERR1", "arm": "primary_amplicon"},
            {
                "run_accession": "ERR1",
                "output_path": str(reads),
                "output_sha256": sha256_file(reads),
            },
            tmp_path / "results",
            {"threads": 2, "timeout": 10},
        )


def test_failed_run_record_resume_and_tamper(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import json

    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    prep = {"run_accession": "ERR1", "output_path": str(reads), "output_sha256": sha256_file(reads)}
    run = {"run_accession": "ERR1", "arm": "primary_amplicon"}
    settings = {"threads": 2, "timeout": 10, "model": "model"}

    def failed(commands, *, timeout, stdout):
        root = Path(json.loads(Path(commands[0][-1]).read_text())["output"])
        (root / "worker.json").write_text(json.dumps({"exit_code": 1, "error": "test failure"}))
        (root / "run_status.json").write_text(
            json.dumps({"schema_version": 1, "status": "execution_failed"})
        )
        raise RuntimeError("test failure")

    monkeypatch.setattr("muc_one_span.clinical_runner.run_tool_pipeline", failed)
    record = run_case(run, prep, tmp_path / "results", settings)
    assert record["status"] == "execution_failed"
    assert record["events"] is None
    assert validate_record(record)
    assert run_case(run, prep, tmp_path / "results", settings, resume=True) == record
    Path(record["outputs"]["worker.json"]["path"]).write_text("{}")
    assert not validate_record(record)
    with pytest.raises(ValueError, match="existing"):
        run_case(run, prep, tmp_path / "results", settings, resume=True)


def test_timeout_without_worker_is_failed_and_not_resumable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")

    def timeout(*args, **kwargs):
        raise TimeoutError("lifecycle budget")

    monkeypatch.setattr("muc_one_span.clinical_runner.run_tool_pipeline", timeout)
    result = run_case(
        {"run_accession": "ERR1", "arm": "primary_amplicon"},
        {"run_accession": "ERR1", "output_path": str(reads), "output_sha256": sha256_file(reads)},
        tmp_path / "out",
        {"threads": 2, "timeout": 1, "model": "unused"},
    )
    assert result["status"] == "execution_failed"
    assert not validate_record(result)
    assert result["exit_code"] is None
    assert "TimeoutError" in result["error"]


def test_wrong_preparation_identity_is_rejected(tmp_path: Path):
    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="accession mismatch"):
        run_case(
            {"run_accession": "ERR1", "arm": "primary_amplicon"},
            {
                "run_accession": "ERR2",
                "output_path": str(reads),
                "output_sha256": sha256_file(reads),
            },
            tmp_path / "out",
            {"threads": 2, "timeout": 10, "model": "unused"},
        )


def test_baseline_verification_rejects_source_drift(tmp_path: Path):
    from muc_one_span.clinical_provenance import verify_environment

    root = tmp_path / "src" / "muc_one_span"
    root.mkdir(parents=True)
    (root / "pipeline.py").write_text("changed")
    with pytest.raises(ValueError, match="source"):
        verify_environment({"caller_sources": {"pipeline.py": "0" * 64}}, tmp_path)


@pytest.mark.parametrize("failure", ["malformed", "timeout_nocall", "interrupted"])
def test_worker_failure_is_terminal(tmp_path: Path, monkeypatch, failure):
    import json

    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")

    def execute(commands, **kwargs):
        root = Path(json.loads(Path(commands[0][-1]).read_text())["output"])
        if failure == "malformed":
            (root / "worker.json").write_text("{")
        elif failure == "timeout_nocall":
            (root / "worker.json").write_text(json.dumps({"exit_code": 1}))
            (root / "run_status.json").write_text(
                json.dumps({"schema_version": 1, "status": "insufficient_evidence"})
            )
            raise TimeoutError("cleanup deadline")
        else:
            raise KeyboardInterrupt

    monkeypatch.setattr("muc_one_span.clinical_runner.run_tool_pipeline", execute)
    args = (
        {"run_accession": "ERR1", "arm": "primary_amplicon"},
        {"run_accession": "ERR1", "output_path": str(reads), "output_sha256": sha256_file(reads)},
        tmp_path / "out",
        {"threads": 2, "timeout": 10, "model": "unused"},
    )
    if failure == "interrupted":
        with pytest.raises(KeyboardInterrupt):
            run_case(*args)
    else:
        run_case(*args)
    record = json.loads((tmp_path / "out/ERR1/benchmark_record.json").read_text())
    assert record["status"] == "execution_failed"
    assert record["ended_at"]


def test_successful_record_is_bound_to_identity_and_all_evidence(tmp_path: Path, monkeypatch):
    import json
    from copy import deepcopy

    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")

    def execute(commands, **kwargs):
        root = Path(json.loads(Path(commands[0][-1]).read_text())["output"])
        data = {
            "classifications": {"allele_1": {"structure": "A", "mutations": []}},
            "alleles": {"allele_1": {"length": 1, "canonical_repeats": 0}},
        }
        for name, payload in {
            "worker.json": {"exit_code": 0},
            "summary.json": data,
            "run_status.json": {"schema_version": 1, "status": "completed"},
            "run_configuration.json": {"input_sha256": sha256_file(reads)},
        }.items():
            (root / name).write_text(json.dumps(payload))
        (root / "consensus_allele_1.fa").write_text(">allele_1\nACGT\n")

    monkeypatch.setattr("muc_one_span.clinical_runner.run_tool_pipeline", execute)
    args = (
        {"run_accession": "ERR1", "arm": "primary_amplicon"},
        {"run_accession": "ERR1", "output_path": str(reads), "output_sha256": sha256_file(reads)},
        tmp_path / "out",
        {"threads": 2, "timeout": 10, "model": "unused"},
    )
    record = run_case(*args)
    assert record["status"] == "completed"
    assert validate_record(record)
    assert run_case(*args, resume=True) == record
    spoof = deepcopy(record)
    spoof["run_accession"] = "ERR2"
    assert not validate_record(spoof)
    missing = deepcopy(record)
    del missing["outputs"]["consensus_allele_1.fa"]
    assert not validate_record(missing)


def test_environment_verifies_import_location_and_cpu_setting(tmp_path: Path, monkeypatch):
    import sys

    from muc_one_span import clinical_provenance
    from muc_one_span.clinical_provenance import installed_packages, verify_environment

    environment = {
        "caller_sources": {},
        "python": sys.version,
        "packages": installed_packages(),
        "execution_environment": {"CUDA_VISIBLE_DEVICES": ""},
    }
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    with pytest.raises(ValueError, match="imported caller"):
        verify_environment(environment, tmp_path)
    checkout = Path(clinical_provenance.__file__).resolve().parents[2]
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    with pytest.raises(ValueError, match="execution environment"):
        verify_environment(environment, checkout)


def test_environment_rejects_installed_package_drift(monkeypatch):
    import sys

    from muc_one_span import clinical_provenance
    from muc_one_span.clinical_provenance import verify_environment

    checkout = Path(clinical_provenance.__file__).resolve().parents[2]
    environment = {
        "caller_sources": {},
        "python": sys.version,
        "packages": {"invented": "0"},
        "execution_environment": {"CUDA_VISIBLE_DEVICES": ""},
    }
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    with pytest.raises(ValueError, match="installed packages"):
        verify_environment(environment, checkout)


def test_engine_and_assay_come_from_hashed_settings(tmp_path: Path, monkeypatch):
    import contextlib
    import json

    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    seen: list[list[str]] = []

    def execute(commands, **kwargs):
        seen.append(json.loads(Path(commands[0][-1]).read_text())["argv"])
        raise RuntimeError("stop after recording argv")

    monkeypatch.setattr("muc_one_span.clinical_runner.run_tool_pipeline", execute)
    run = {"run_accession": "ERR1", "arm": "primary_amplicon"}
    prep = {"run_accession": "ERR1", "output_path": str(reads), "output_sha256": sha256_file(reads)}
    base = {"threads": 2, "timeout": 10, "model": "unused"}
    for root, extra in (("ladder", {}), ("hybrid", {"engine": "hybrid", "assay": "genomic"})):
        with contextlib.suppress(RuntimeError):
            run_case(run, prep, tmp_path / root, {**base, **extra})
    assert "--engine" not in seen[0]
    assert seen[1][-4:] == ["--engine", "hybrid", "--assay", "genomic"]
