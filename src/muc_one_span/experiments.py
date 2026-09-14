"""Strict, inventoried MucOneUp amplicon experiments; no caller inference."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import time
from pathlib import Path
from typing import Any

from muc_one_span.config import _bundled_repeats_path, load_repeat_dictionary
from muc_one_span.evaluation.artifacts import discover_input
from muc_one_span.evaluation.truth import load_truth
from muc_one_span.tools import _clean_path_for_externals, run_tool, run_tool_iter

SIMULATOR_PLATFORMS = {"hifi": "pacbio", "ont": "ont"}


def _object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _invalid_constant(value: str) -> None:
    raise ValueError(f"nonfinite JSON value: {value}")


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(
        path.read_text(), object_pairs_hook=_object, parse_constant=_invalid_constant
    )
    if not isinstance(value, dict):
        raise ValueError("expected a JSON object")
    return value


def _keys(value: dict[str, Any], allowed: set[str]) -> None:
    if unknown := value.keys() - allowed:
        raise ValueError(f"unknown fields: {sorted(unknown)}")


def _integer(value: Any, minimum: int, label: str) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _path(value: Any, base: Path) -> str:
    if not isinstance(value, str) or not value or "\x00" in value:
        raise ValueError("paths must be nonempty strings")
    return str((base / value).resolve())


def _case(value: Any) -> dict[str, Any]:
    fields = {"sample", "platform", "seed", "lengths", "mutation", "targets", "requested_templates"}
    if not isinstance(value, dict) or set(value) != fields:
        raise ValueError(f"each case requires exactly {sorted(fields)}")
    if not isinstance(value["sample"], str) or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_.-]*", value["sample"]
    ):
        raise ValueError("sample must be a safe path component without spaces")
    if value["platform"] not in ("hifi", "ont"):
        raise ValueError("platform must be hifi or ont")
    _integer(value["seed"], 0, "seed")
    _integer(value["requested_templates"], 1, "requested_templates")
    lengths = value["lengths"]
    if not isinstance(lengths, list) or len(lengths) != 2:
        raise ValueError("lengths must contain two positive integers")
    for length in lengths:
        _integer(length, 1, "length")
    mutation, targets = value["mutation"], value["targets"]
    if mutation is not None and (
        not isinstance(mutation, str) or not re.fullmatch(r"[A-Za-z0-9_]+", mutation)
    ):
        raise ValueError("mutation must be null or a single mutation name")
    if not isinstance(targets, list) or bool(targets) != (mutation is not None):
        raise ValueError("mutation requires targets; normal cases use null and []")
    seen = set()
    for target in targets:
        if not isinstance(target, list) or len(target) != 2:
            raise ValueError("targets must contain [haplotype, repeat] pairs")
        hap, repeat = (_integer(v, 1, "target") for v in target)
        if hap > 2 or repeat > lengths[hap - 1] or (hap, repeat) in seen:
            raise ValueError("duplicate or out-of-range mutation target")
        seen.add((hap, repeat))
    return value


def load_design(path: Path) -> dict[str, Any]:
    """Load strict schema 1; resolve design paths relative to the JSON file."""
    value = _read_json(path)
    _keys(value, {"schema_version", "cases", "platform_configs", "repeat_dictionary"})
    if type(value.get("schema_version")) is not int or value["schema_version"] != 1:
        raise ValueError("schema_version must be integer 1")
    if not isinstance(value.get("cases"), list) or not value["cases"]:
        raise ValueError("cases must be a nonempty array")
    cases = [_case(c) for c in value["cases"]]
    if len({c["sample"] for c in cases}) != len(cases):
        raise ValueError("duplicate sample names")
    configs = value.get("platform_configs", {})
    if not isinstance(configs, dict):
        raise ValueError("platform_configs must be an object")
    _keys(configs, set(SIMULATOR_PLATFORMS))
    value["platform_configs"] = {k: _path(v, path.resolve().parent) for k, v in configs.items()}
    if "repeat_dictionary" in value:
        value["repeat_dictionary"] = _path(value["repeat_dictionary"], path.resolve().parent)
    return value


def case_commands(
    case: dict[str, Any], config: Path, output: Path, executable: str
) -> list[list[str]]:
    """Build exact simulator argument vectors with one diploid FASTA output."""
    name = case["sample"]
    base = [executable, "--config", str(config)]
    simulate = [
        *base,
        "simulate",
        "--out-base",
        name,
        "--out-dir",
        str(output),
        "--num-haplotypes",
        "2",
        "--fixed-lengths",
        str(case["lengths"][0]),
        "--fixed-lengths",
        str(case["lengths"][1]),
        "--output-structure",
        "--seed",
        str(case["seed"]),
    ]
    if case["mutation"] is not None:
        simulate += ["--mutation-name", case["mutation"]]
        for hap, repeat in case["targets"]:
            simulate += ["--mutation-targets", f"{hap},{repeat}"]
    reads = [
        *base,
        "reads",
        "amplicon",
        str(output / f"{name}.001.simulated.fa"),
        "--out-dir",
        str(output),
        "--out-base",
        f"{name}_reads",
        "--coverage",
        str(case["requested_templates"]),
        "--seed",
        str(case["seed"]),
        "--platform",
        SIMULATOR_PLATFORMS[case["platform"]],
    ]
    return [simulate, reads]


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_artifacts(hashes: dict[str, str]) -> None:
    for artifact, digest in hashes.items():
        if _hash(Path(artifact)) != digest:
            raise ValueError(f"artifact changed during execution: {artifact}")


def _write(path: Path, data: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def count_usable_records(path: Path) -> int:
    """Count primary BAM records or valid FASTQ records, preserving name collisions."""
    count = 0
    if path.suffix == ".bam":
        run_tool(["samtools", "quickcheck", str(path)])
        for line in run_tool_iter(["samtools", "view", "-F", "2304", str(path)]):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("malformed SAM record")
            sequence, quality = fields[9:11]
            if sequence != "*" and quality != "*" and len(sequence) == len(quality):
                count += 1
    else:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            while header := handle.readline():
                sequence, plus, quality = (handle.readline().rstrip("\r\n") for _ in range(3))
                if (
                    not header.startswith("@")
                    or not plus.startswith("+")
                    or not sequence
                    or len(sequence) != len(quality)
                ):
                    raise ValueError("malformed FASTQ record")
                if any(ord(c) < 33 or ord(c) > 126 for c in quality):
                    raise ValueError("invalid FASTQ quality")
                count += 1
    if not count:
        raise ValueError("no usable primary records with sequence and quality")
    return count


def _config_evidence(path: Path, platform: str) -> dict[str, Any]:
    config = _read_json(path)
    section = "pacbio_params" if platform == "hifi" else "ont_amplicon_params"
    parameters = config.get(section, {})
    if not isinstance(parameters, dict):
        raise ValueError(f"{section} must be an object")
    model = parameters.get("model_file")
    if not isinstance(model, str) or not model:
        raise ValueError(f"{section}.model_file must be explicit for reproducible experiments")
    model_path = Path(_path(model, path.parent))
    files: dict[str, str] = {str(path): _hash(path), str(model_path): _hash(model_path)}

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            for nested in value.values():
                visit(nested)
        elif isinstance(value, list):
            for nested in value:
                visit(nested)
        elif isinstance(value, str) and value and "\x00" not in value:
            candidate = path.parent / value
            try:
                if candidate.is_file():
                    files[str(candidate.resolve())] = _hash(candidate)
            except OSError:
                pass  # Sequence strings and tool command strings are not filesystem paths.

    visit(config)
    return {
        "path": str(path),
        "sha256": files[str(path)],
        "model": str(model_path),
        "artifact_hashes": files,
        "snapshot": config,
        "hash_scope": "config, explicit model and existing file-valued config entries; no environment lock",
    }


def run_experiment(
    design_path: Path,
    output_dir: Path,
    *,
    config: Path | None = None,
    executable: str = "muconeup",
    dry_run: bool = False,
    fallback_config: Path | None = None,
) -> dict[str, Any]:
    """Execute every case, preserving failures and refusing existing case results.

    Explicit config wins over platform configs, then fallback_config. Dry runs
    only validate the design and return planned commands, without tool execution.
    """
    design = load_design(design_path)
    if "/" in executable:
        executable = str(Path(executable).resolve())
    found_executable = shutil.which(
        executable, path=_clean_path_for_externals(os.environ.get("PATH", ""))
    )
    resolved_executable = str(Path(found_executable).absolute()) if found_executable else None
    if resolved_executable is not None:
        executable = resolved_executable
    output_dir = output_dir.resolve()
    selected = []
    for case in design["cases"]:
        cfg = config or design["platform_configs"].get(case["platform"]) or fallback_config
        if cfg is None:
            raise ValueError(f"no MucOneUp config for {case['platform']}")
        selected.append(Path(cfg).resolve())
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "design": design,
        "design_path": str(design_path.resolve()),
        "design_sha256": _hash(design_path),
        "requested_templates_are_not_retained_records": True,
        "read_source_truth": "unavailable: MucOneUp amplicon source tracking unsupported",
        "cases": [
            {
                **case,
                "status": "planned",
                "commands": [],
                "planned_commands": case_commands(
                    case, cfg, output_dir / case["sample"], executable
                ),
            }
            for case, cfg in zip(design["cases"], selected, strict=True)
        ],
    }
    if dry_run:
        return manifest
    if (output_dir / "generation_manifest.json").exists() or any(
        (output_dir / c["sample"]).exists() for c in design["cases"]
    ):
        raise FileExistsError(
            "existing experiment manifest or case directory; select a new output directory"
        )
    if (output_dir / "inventory.json").exists():
        raise FileExistsError("existing experiment inventory")
    output_dir.mkdir(parents=True, exist_ok=True)
    inventory = [
        {
            "sample": c["sample"],
            "platform": c["platform"],
            "input": str(output_dir / c["sample"] / "__generation_failed__.fastq"),
            "truth_dir": str(output_dir / c["sample"]),
            "generation_status": "planned",
        }
        for c in design["cases"]
    ]
    manifest_path = output_dir / "generation_manifest.json"

    def save() -> None:
        _write(manifest_path, manifest)
        _write(output_dir / "inventory.json", inventory)

    save()
    version_cmd = [executable, "--version"]
    manifest["version_command"] = version_cmd
    manifest["version_result"] = {"status": "running"}
    save()
    try:
        manifest["simulator_version"] = run_tool(version_cmd).strip()
        manifest["version_result"] = {
            "status": "completed",
            "stdout": manifest["simulator_version"],
        }
        manifest["simulator_executable"] = resolved_executable
        manifest["simulator_executable_sha256"] = (
            _hash(Path(resolved_executable)) if resolved_executable else None
        )
        dictionary_path = Path(design.get("repeat_dictionary", _bundled_repeats_path()))
        rd = load_repeat_dictionary(dictionary_path)
        manifest["repeat_dictionary"] = {
            "path": str(dictionary_path),
            "sha256": _hash(dictionary_path),
        }
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        manifest["setup_error"] = str(exc)
        if manifest["version_result"]["status"] == "running":
            manifest["version_result"] = {"status": "execution_failed", "error": str(exc)}
        for row, entry in zip(manifest["cases"], inventory, strict=True):
            row.update(status="execution_failed", error=str(exc))
            entry["generation_status"] = "execution_failed"
        save()
        return manifest
    manifest["version_command"] = version_cmd
    save()
    for row, entry, cfg in zip(manifest["cases"], inventory, selected, strict=True):
        destination = output_dir / row["sample"]
        destination.mkdir()
        started = time.monotonic()
        try:
            evidence = _config_evidence(cfg, row["platform"])
            row["config"] = evidence
            _write(destination / "simulator_config.snapshot.json", evidence["snapshot"])
            if row["mutation"] is not None and row["mutation"] not in rd.mutations:
                raise ValueError(f"mutation absent from truth dictionary: {row['mutation']}")
            for command in row["planned_commands"]:
                _verify_artifacts(evidence["artifact_hashes"])
                result: dict[str, Any] = {
                    "argv": command,
                    "cwd": str(cfg.parent),
                    "status": "running",
                }
                row["commands"].append(result)
                save()
                try:
                    result["stdout"] = run_tool(command, cwd=str(cfg.parent))
                    result["status"] = "completed"
                except (OSError, RuntimeError) as exc:
                    result.update(status="execution_failed", error=str(exc))
                    raise
                _verify_artifacts(evidence["artifact_hashes"])
            truth = load_truth(destination, rd)
            if (
                len(truth.haplotypes) != 2
                or [len(h.structure) for h in truth.haplotypes] != row["lengths"]
            ):
                raise ValueError("generated truth disagrees with requested diploid lengths")
            actual_events = sorted(
                (i, event.repeat_index, event.name)
                for i, hap in enumerate(truth.haplotypes, 1)
                for event in hap.events
            )
            expected_events = sorted((h, r, row["mutation"]) for h, r in row["targets"])
            if actual_events != expected_events:
                raise ValueError("generated truth disagrees with requested mutation targets")
            input_path = discover_input(destination)
            row["usable_records"] = count_usable_records(input_path)
            row.update(
                status="completed",
                input=str(input_path),
                input_sha256=_hash(input_path),
                truth_hashes=truth.hashes,
                truth_warnings=list(truth.warnings),
            )
            entry["input"] = str(input_path)
        except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
            row.update(status="execution_failed", error=str(exc))
        row["wall_seconds"] = time.monotonic() - started
        row["output_hashes"] = {
            str(p.relative_to(destination)): _hash(p)
            for p in sorted(destination.rglob("*"))
            if p.is_file()
        }
        entry["generation_status"] = row["status"]
        save()
    return manifest
