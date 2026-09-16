"""Auditable ENA inventory retrieval and identity-preserving FASTQ preparation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
import shutil
import tempfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

SCHEMA_VERSION = 1
ENA_FIELDS = (
    "run_accession,sample_accession,experiment_accession,sample_alias,scientific_name,"
    "instrument_platform,instrument_model,library_strategy,library_source,library_selection,"
    "library_layout,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5"
)
REQUIRED_ENA_COLUMNS = {
    "run_accession",
    "sample_accession",
    "experiment_accession",
    "sample_alias",
    "instrument_platform",
    "library_strategy",
    "fastq_ftp",
    "fastq_bytes",
    "fastq_md5",
}


def fetch_inventory(study: str) -> dict[str, Any]:
    """Fetch and convert the complete ENA read-run report for ``study``."""
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={study}&result=read_run&fields={ENA_FIELDS}&format=tsv&download=false"
    )
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            text = response.read().decode("utf-8")
    except (OSError, UnicodeError, urllib.error.URLError) as exc:
        raise ValueError(f"Could not retrieve ENA inventory for {study}: {exc}") from exc
    retrieved = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    return inventory_from_tsv(text, study=study, source_url=url, retrieved_at=retrieved)


def inventory_from_tsv(
    text: str, *, study: str, source_url: str, retrieved_at: str
) -> dict[str, Any]:
    """Convert a saved ENA TSV report into the versioned inventory schema."""
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    columns = set(reader.fieldnames or ())
    missing = REQUIRED_ENA_COLUMNS - columns
    if missing:
        raise ValueError(f"missing ENA columns: {', '.join(sorted(missing))}")
    runs: list[dict[str, Any]] = []
    try:
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("ENA row has the wrong number of columns")
            runs.append(_convert_row(dict(row)))
    except csv.Error as exc:
        raise ValueError(f"malformed ENA metadata: {exc}") from exc
    data = {
        "schema_version": SCHEMA_VERSION,
        "study_accession": study,
        "retrieved_at": retrieved_at,
        "source_url": source_url,
        "runs": runs,
    }
    return validate_inventory(data)


def _convert_row(row: dict[str, str]) -> dict[str, Any]:
    alias = row["sample_alias"].strip()
    strategy = row["library_strategy"].strip().upper()
    target = alias.rsplit("_", 1)[-1].upper()
    if strategy == "AMPLICON" and target == "MUC1":
        arm, reason = "primary_amplicon", "MUC1-targeted PCR amplicon"
    elif strategy == "WGS" and target == "MUC1":
        arm, reason = "secondary_wgs", "MUC1-locus reads submitted from WGS"
    else:
        arm, reason = "excluded", f"off-target or unsupported library ({target}, {strategy})"
    files = _file_entries(row["fastq_ftp"], row["fastq_bytes"], row["fastq_md5"])
    biological = alias.split("_PCR_", 1)[0].split("_WGS_", 1)[0]
    return {
        "run_accession": row["run_accession"].strip(),
        "sample_accession": row["sample_accession"].strip(),
        "experiment_accession": row["experiment_accession"].strip(),
        "alias": alias,
        "biological_sample": biological,
        "biological_sample_basis": "provisional alias prefix; confirm in truth ledger",
        "identity_resolution": "alias_only",
        "arm": arm,
        "inclusion_reason": reason,
        "platform": row["instrument_platform"].strip(),
        "files": files,
        "ena": row,
    }


def _file_entries(urls: str, sizes: str, md5s: str) -> list[dict[str, Any]]:
    url_items, size_items, md5_items = urls.split(";"), sizes.split(";"), md5s.split(";")
    if not urls or not (len(url_items) == len(size_items) == len(md5_items)):
        raise ValueError("ENA FASTQ URL, byte, and MD5 fields do not align")
    entries = []
    for url, size, md5 in zip(url_items, size_items, md5_items, strict=True):
        normalized = url if "://" in url else f"https://{url}"
        try:
            byte_count = int(size)
        except ValueError as exc:
            raise ValueError(f"invalid ENA byte count: {size}") from exc
        if byte_count < 0 or len(md5) != 32:
            raise ValueError("invalid ENA file size or MD5")
        entries.append({"url": normalized, "bytes": byte_count, "md5": md5.lower()})
    return entries


def validate_inventory(data: dict[str, Any]) -> dict[str, Any]:
    """Validate inventory structure while permitting replicated biological identities."""
    if data.get("schema_version") != SCHEMA_VERSION or not isinstance(data.get("runs"), list):
        raise ValueError("invalid clinical inventory schema")
    for field in ("study_accession", "retrieved_at", "source_url"):
        if not isinstance(data.get(field), str) or not data[field]:
            raise ValueError(f"missing inventory {field}")
    if not data["runs"]:
        raise ValueError("empty ENA inventory")
    seen: set[str] = set()
    for run in data["runs"]:
        if not isinstance(run, dict):
            raise ValueError("run must be an object")
        accession = _accession(run.get("run_accession"))
        for field in (
            "sample_accession",
            "experiment_accession",
            "alias",
            "biological_sample",
            "inclusion_reason",
            "platform",
        ):
            if not isinstance(run.get(field), str) or not run[field]:
                raise ValueError(f"missing run {field}")
        if accession in seen:
            raise ValueError(f"duplicate run_accession: {accession}")
        seen.add(accession)
        if run.get("arm") not in {"primary_amplicon", "secondary_wgs", "excluded"}:
            raise ValueError(f"invalid arm for {accession}")
        if not isinstance(run.get("ena"), dict) or not isinstance(run.get("files"), list):
            raise ValueError(f"missing preserved metadata or files for {accession}")
        _validate_files(run["files"])
    return data


def _accession(value: Any) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[EDS]RR[0-9]+", value) is None:
        raise ValueError("invalid run accession")
    return value


def _validate_files(files: Any) -> None:
    if not isinstance(files, list) or not files:
        raise ValueError("run requires nonempty files")
    names: set[str] = set()
    for item in files:
        if not isinstance(item, dict):
            raise ValueError("file entry must be an object")
        url = urlparse(str(item.get("url", "")))
        name = Path(url.path).name
        if url.scheme != "https" or not url.netloc or not name or name in names:
            raise ValueError("invalid or colliding HTTPS file URL")
        names.add(name)
        size = item.get("bytes")
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            raise ValueError("invalid file byte count")
        if re.fullmatch(r"[a-fA-F0-9]{32}", str(item.get("md5", ""))) is None:
            raise ValueError("invalid file MD5")


def prepare_run(run: dict[str, Any], data_root: Path) -> dict[str, Any]:
    """Download, validate, and decompress one run without changing read content."""
    accession = _accession(run["run_accession"])
    _validate_files(run.get("files"))
    run_dir = data_root / accession
    run_dir.mkdir(parents=True, exist_ok=True)
    prepared_files = []
    inputs: list[Path] = []
    for file_data in run["files"]:
        name = Path(urlparse(str(file_data["url"])).path).name
        if not name or name in {".", ".."}:
            raise ValueError("download URL has no safe filename")
        destination = run_dir / name
        _ensure_download(file_data, destination)
        prepared_files.append(
            {**file_data, "sha256": _hash_file(destination, "sha256"), "path": str(destination)}
        )
        inputs.append(destination)
    output = run_dir / f"{accession}.fastq"
    while output in inputs:
        output = output.with_name(output.stem + ".prepared.fastq")
    with tempfile.NamedTemporaryFile(dir=run_dir, suffix=".preparing", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with temporary.open("wb") as target:
            for source in inputs:
                opener = gzip.open if source.suffix == ".gz" else Path.open
                with opener(source, "rb") as stream:
                    shutil.copyfileobj(stream, target)
        stats = _inspect_fastq(temporary)
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    combined_input_hash = hashlib.sha256()
    for source in inputs:
        with source.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                combined_input_hash.update(chunk)
    return {
        "schema_version": SCHEMA_VERSION,
        "run_accession": accession,
        "input_sha256": combined_input_hash.hexdigest(),
        "output_sha256": _hash_file(output, "sha256"),
        **stats,
        "excluded_reads": 0,
        "preprocessing": "validated identity-preserving decompression",
        "output_path": str(output),
        "files": prepared_files,
    }


def _ensure_download(file_data: dict[str, Any], destination: Path) -> None:
    expected_size, expected_md5 = int(file_data["bytes"]), str(file_data["md5"]).lower()
    if destination.exists():
        _verify_file(destination, expected_size, expected_md5)
        return
    partial = destination.with_name(destination.name + ".part")
    offset = partial.stat().st_size if partial.exists() else 0
    if offset == expected_size:
        try:
            _verify_file(partial, expected_size, expected_md5)
        except ValueError:
            partial.unlink()
            offset = 0
        else:
            partial.replace(destination)
            return
    if offset > expected_size:
        partial.unlink()
        offset = 0
    request = urllib.request.Request(str(file_data["url"]))
    if offset:
        request.add_header("Range", f"bytes={offset}-")
    try:
        started = time.monotonic()
        with urllib.request.urlopen(request, timeout=60) as response:
            mode = "ab" if offset and getattr(response, "status", 206) == 206 else "wb"
            written = offset if mode == "ab" else 0
            with partial.open(mode) as target:
                while chunk := response.read(1024 * 1024):
                    if time.monotonic() - started > 1800:
                        raise ValueError("download exceeded 1800-second transfer budget")
                    written += len(chunk)
                    if written > expected_size:
                        raise ValueError("download exceeds expected byte count")
                    target.write(chunk)
    except (OSError, urllib.error.URLError) as exc:
        raise ValueError(f"download failed for {file_data['url']}: {exc}") from exc
    _verify_file(partial, expected_size, expected_md5)
    partial.replace(destination)


def _verify_file(path: Path, expected_size: int, expected_md5: str) -> None:
    if path.stat().st_size != expected_size:
        raise ValueError(f"byte count mismatch for {path}")
    if _hash_file(path, "md5") != expected_md5:
        raise ValueError(f"MD5 mismatch for {path}")


def _hash_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _inspect_fastq(path: Path) -> dict[str, Any]:
    lengths: list[int] = []
    first_header = last_header = ""
    read_ids: set[str] = set()
    id_digest = hashlib.sha256()
    with path.open("rb") as stream:
        while header := stream.readline():
            sequence, plus, quality = stream.readline(), stream.readline(), stream.readline()
            if (
                not (sequence and plus and quality)
                or not header.startswith(b"@")
                or not plus.startswith(b"+")
            ):
                raise ValueError("invalid FASTQ record structure")
            sequence, quality = sequence.rstrip(b"\r\n"), quality.rstrip(b"\r\n")
            if len(sequence) != len(quality):
                raise ValueError("FASTQ sequence and quality length differ")
            header_text = header.rstrip(b"\r\n").decode("utf-8")
            parts = header_text[1:].split(maxsplit=1)
            read_id = parts[0] if parts else ""
            if not read_id or read_id in read_ids:
                raise ValueError(f"duplicate or empty FASTQ read ID: {read_id}")
            read_ids.add(read_id)
            id_digest.update(read_id.encode() + b"\n")
            if not first_header:
                first_header = header_text
            last_header = header_text
            lengths.append(len(sequence))
    if not lengths:
        raise ValueError("FASTQ contains no reads")
    ordered = sorted(lengths)

    def quantile(percent: int) -> int:
        return ordered[(len(ordered) - 1) * percent // 100]

    return {
        "input_reads": len(lengths),
        "output_reads": len(lengths),
        "read_lengths": {
            "min": ordered[0],
            "p25": quantile(25),
            "median": quantile(50),
            "p75": quantile(75),
            "max": ordered[-1],
            "total": sum(ordered),
        },
        "headers": {"first": first_header, "last": last_header},
        "read_id_sha256": id_digest.hexdigest(),
    }
