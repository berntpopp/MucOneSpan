"""Tests for auditable ENA inventory conversion and read preparation."""

from __future__ import annotations

import gzip
import hashlib
import io
from pathlib import Path

import pytest

from muc_one_span.clinical_data import inventory_from_tsv, prepare_run

HEADER = "\t".join(
    (
        "run_accession",
        "sample_accession",
        "experiment_accession",
        "sample_alias",
        "instrument_platform",
        "library_strategy",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
    )
)


def _row(run: str = "ERR1", alias: str = "HG001_PCR_MUC1", payload: bytes = b"x") -> str:
    strategy = "WGS" if "_WGS_" in alias else "AMPLICON"
    return "\t".join(
        (
            run,
            "SAMEA1",
            "ERX1",
            alias,
            "OXFORD_NANOPORE",
            strategy,
            f"ftp.example/{run}.fastq.gz",
            str(len(payload)),
            hashlib.md5(payload).hexdigest(),
        )
    )


def _inventory(*rows: str) -> dict:
    return inventory_from_tsv(
        HEADER + "\n" + "\n".join(rows) + "\n",
        study="PRJEB92208",
        source_url="https://example/api",
        retrieved_at="2026-09-15T00:00:00Z",
    )


def test_offline_conversion_preserves_metadata_ids_and_repeated_identity() -> None:
    data = _inventory(_row("ERR1", "HG001_PCR_MUC1"), _row("ERR2", "HG001_WGS_MUC1"))
    first, second = data["runs"]
    assert first["ena"]["sample_alias"] == "HG001_PCR_MUC1"
    assert first["sample_accession"] == "SAMEA1"
    assert first["biological_sample"] == second["biological_sample"] == "HG001"
    assert [first["arm"], second["arm"]] == ["primary_amplicon", "secondary_wgs"]
    assert first["files"] == [
        {
            "url": "https://ftp.example/ERR1.fastq.gz",
            "bytes": 1,
            "md5": "9dd4e461268c8034f5c8564e155c67a6",
        }
    ]


def test_offline_conversion_classifies_off_target_runs_without_dropping_them() -> None:
    run = _inventory(_row(alias="HG001_PCR_ACAN"))["runs"][0]
    assert run["arm"] == "excluded"
    assert "ACAN" in run["inclusion_reason"]


@pytest.mark.parametrize(
    "text,match",
    [
        (HEADER + "\nERR1\ttoo-few\n", "columns"),
        (HEADER.replace("fastq_md5", "other") + "\n" + _row(), "missing ENA columns"),
    ],
)
def test_offline_conversion_rejects_malformed_metadata(text: str, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        inventory_from_tsv(text, study="P", source_url="u", retrieved_at="t")


def test_validate_inventory_rejects_duplicate_run_accessions() -> None:
    with pytest.raises(ValueError, match="duplicate run_accession"):
        _inventory(_row(), _row())


def _prepared_run(payload: bytes, *, size: int | None = None, md5: str | None = None) -> dict:
    return {
        "run_accession": "ERR1",
        "files": [
            {
                "url": "https://example/reads.fastq.gz",
                "bytes": len(payload) if size is None else size,
                "md5": hashlib.md5(payload).hexdigest() if md5 is None else md5,
            }
        ],
    }


def test_prepare_run_replays_local_gzip_and_preserves_reads(tmp_path: Path) -> None:
    raw = b"@r1 comment\nACGT\n+\n!!!!\n@r2\nAA\n+\n##\n"
    compressed = gzip.compress(raw, mtime=0)
    source = tmp_path / "ERR1" / "reads.fastq.gz"
    source.parent.mkdir()
    source.write_bytes(compressed)
    record = prepare_run(_prepared_run(compressed), tmp_path)
    assert (tmp_path / "ERR1" / "ERR1.fastq").read_bytes() == raw
    assert record["input_reads"] == record["output_reads"] == 2
    assert record["excluded_reads"] == 0
    assert record["read_lengths"] == {
        "min": 2,
        "p25": 2,
        "median": 2,
        "p75": 2,
        "max": 4,
        "total": 6,
    }
    assert record["headers"] == {"first": "@r1 comment", "last": "@r2"}
    assert record["input_sha256"] == hashlib.sha256(compressed).hexdigest()
    assert record["output_sha256"] == hashlib.sha256(raw).hexdigest()


@pytest.mark.parametrize("size,md5,match", [(999, None, "byte count"), (None, "0" * 32, "MD5")])
def test_prepare_run_rejects_stale_local_file(
    tmp_path: Path, size: int | None, md5: str | None, match: str
) -> None:
    payload = gzip.compress(b"@r\nA\n+\n!\n", mtime=0)
    path = tmp_path / "ERR1" / "reads.fastq.gz"
    path.parent.mkdir()
    path.write_bytes(payload)
    with pytest.raises(ValueError, match=match):
        prepare_run(_prepared_run(payload, size=size, md5=md5), tmp_path)


def test_prepare_run_rejects_invalid_fastq(tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nAC\n+\n!\n", mtime=0)
    path = tmp_path / "ERR1" / "reads.fastq.gz"
    path.parent.mkdir()
    path.write_bytes(payload)
    with pytest.raises(ValueError, match="quality length"):
        prepare_run(_prepared_run(payload), tmp_path)


def test_prepare_run_rejects_duplicate_read_ids(tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nA\n+\n!\n@r comment\nC\n+\n#\n", mtime=0)
    path = tmp_path / "ERR1" / "reads.fastq.gz"
    path.parent.mkdir()
    path.write_bytes(payload)
    with pytest.raises(ValueError, match="duplicate"):
        prepare_run(_prepared_run(payload), tmp_path)


class _Response(io.BytesIO):
    status = 206

    def __enter__(self) -> _Response:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


def test_prepare_run_resumes_part_download(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nA\n+\n!\n", mtime=0)
    partial = tmp_path / "ERR1" / "reads.fastq.gz.part"
    partial.parent.mkdir()
    partial.write_bytes(payload[:5])
    seen_range: list[str | None] = []

    def open_remainder(request: object, timeout: int) -> _Response:
        seen_range.append(request.get_header("Range"))  # type: ignore[attr-defined]
        assert timeout == 60
        return _Response(payload[5:])

    monkeypatch.setattr("muc_one_span.clinical_data.urllib.request.urlopen", open_remainder)
    record = prepare_run(_prepared_run(payload), tmp_path)
    assert seen_range == ["bytes=5-"]
    assert record["output_reads"] == 1


def test_prepare_run_keeps_truncated_part_for_restart(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    payload = gzip.compress(b"@r\nA\n+\n!\n", mtime=0)
    monkeypatch.setattr(
        "muc_one_span.clinical_data.urllib.request.urlopen",
        lambda request, timeout: _Response(payload[:5]),
    )
    with pytest.raises(ValueError, match="byte count"):
        prepare_run(_prepared_run(payload), tmp_path)
    assert (tmp_path / "ERR1" / "reads.fastq.gz.part").read_bytes() == payload[:5]


def test_plain_source_collision_does_not_destroy_input(tmp_path: Path) -> None:
    raw = b"@r\nAC\n+\n!!\n"
    run = _prepared_run(raw)
    run["files"][0]["url"] = "https://example/ERR1.fastq"
    source = tmp_path / "ERR1" / "ERR1.fastq"
    source.parent.mkdir()
    source.write_bytes(raw)
    record = prepare_run(run, tmp_path)
    assert source.read_bytes() == raw
    assert Path(record["output_path"]).read_bytes() == raw


def test_completed_partial_promoted_offline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = gzip.compress(b"@r\nAC\n+\n!!\n")
    partial = tmp_path / "ERR1" / "reads.fastq.gz.part"
    partial.parent.mkdir()
    partial.write_bytes(payload)
    monkeypatch.setattr(
        "muc_one_span.clinical_data.urllib.request.urlopen",
        lambda *a, **k: pytest.fail("complete valid partial needs no network"),
    )
    assert prepare_run(_prepared_run(payload), tmp_path)["output_reads"] == 1


def test_invalid_preparation_preserves_previous_output(tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nAC\n+\n!\n")
    source = tmp_path / "ERR1" / "reads.fastq.gz"
    source.parent.mkdir()
    source.write_bytes(payload)
    output = source.parent / "ERR1.fastq"
    output.write_bytes(b"previous validated output")
    with pytest.raises(ValueError):
        prepare_run(_prepared_run(payload), tmp_path)
    assert output.read_bytes() == b"previous validated output"


def test_rejects_traversal_accession(tmp_path: Path) -> None:
    run = _prepared_run(b"x")
    run["run_accession"] = "../escape"
    with pytest.raises(ValueError, match="accession"):
        prepare_run(run, tmp_path)


@pytest.mark.parametrize(
    "names", [("ERR1.fastq", "ERR1.prepared.fastq"), ("ERR1.fastq.preparing",)]
)
def test_preparation_preserves_all_colliding_source_names(
    tmp_path: Path, names: tuple[str, ...]
) -> None:
    run = {"run_accession": "ERR1", "files": []}
    root = tmp_path / "ERR1"
    root.mkdir()
    originals = {}
    for index, name in enumerate(names):
        payload = f"@r{index}\nAC\n+\n!!\n".encode()
        (root / name).write_bytes(payload)
        originals[name] = payload
        run["files"].append(
            {
                "url": f"https://example/{name}",
                "bytes": len(payload),
                "md5": hashlib.md5(payload).hexdigest(),
            }
        )
    record = prepare_run(run, tmp_path)
    assert record["output_reads"] == len(names)
    assert all((root / name).read_bytes() == value for name, value in originals.items())
