"""Strict adapter for actual MucOneUp FASTA, structure and statistics artifacts."""

from __future__ import annotations

import hashlib
import json
import shlex
from pathlib import Path
from typing import Any

from muc_one_span.config import RepeatDictionary, _apply_mutation

from .models import Event, TruthHaplotype, TruthSample


class TruthValidationError(ValueError):
    """Missing, ambiguous or inconsistent authoritative simulator truth."""


def fasta_records(path: Path) -> dict[str, str]:
    """Read FASTA without case normalization; reject empty or duplicate records."""
    records: dict[str, str] = {}
    name: str | None = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            parts = line[1:].split()
            if not parts or parts[0] in records:
                raise ValueError(f"{path}: empty or duplicate FASTA identifier")
            name = parts[0]
            records[name] = ""
        elif name is None:
            raise ValueError(f"{path}: sequence before FASTA header")
        else:
            if any(c.isspace() for c in line):
                raise ValueError(f"{path}: whitespace in sequence")
            records[name] += line
    if not records or any(not seq for seq in records.values()):
        raise ValueError(f"{path}: missing or empty FASTA sequence")
    return records


def _one(root: Path, suffix: str) -> Path:
    files = sorted(root.glob(f"*{suffix}"))
    if len(files) != 1:
        raise ValueError(f"expected exactly one {suffix}, found {len(files)}")
    return files[0]


def _option_integer(command: str, flag: str) -> int | None:
    try:
        parts = shlex.split(command)
        value = parts[parts.index(flag) + 1] if flag in parts else None
        return int(value) if value is not None else None
    except (ValueError, IndexError):
        return None


def _read_provenance(root: Path, provenance: dict[str, Any], paths: list[Path]) -> None:
    metadata_paths = sorted(root.glob("*_metadata.tsv"))
    metadata = {}
    if len(metadata_paths) == 1:
        path = metadata_paths[0]
        paths.append(path)
        metadata = dict(
            line.split("\t", 1) for line in path.read_text().splitlines() if "\t" in line
        )
    provenance["read_metadata"] = metadata
    provenance["read_metadata_status"] = (
        "available" if len(metadata_paths) == 1 else "missing_or_ambiguous"
    )
    command = metadata.get("Command", "")
    provenance["nominal_molecules"] = _option_integer(command, "--coverage")
    provenance["retained_records"] = None
    provenance["read_simulation_seed"] = _option_integer(command, "--seed")
    provenance["truth_command_seed"] = _option_integer(provenance.get("command_line", ""), "--seed")
    provenance["platform"] = metadata.get("Read_simulation_technology")
    provenance["model_hashes"] = None


def _load(root: Path, rd: RepeatDictionary) -> TruthSample:
    paths = [
        _one(root, suffix)
        for suffix in (".simulated.fa", ".vntr_structure.txt", ".simulation_stats.json")
    ]
    records = fasta_records(paths[0])
    structures: dict[str, list[str]] = {}
    for line in paths[1].read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        name, chain = line.split("\t")
        if name in structures or not chain:
            raise ValueError("duplicate haplotype or empty structure")
        structures[name] = chain.split("-")
    stats = json.loads(paths[2].read_text())
    hap_stats = stats["haplotype_statistics"]
    expected_ids = [f"haplotype_{i + 1}" for i in range(len(hap_stats))]
    if not expected_ids or set(records) != set(expected_ids) or set(structures) != set(records):
        raise ValueError("FASTA, structure and statistics haplotypes disagree")
    mutation_info = stats.get("mutation_info") or {}
    mutation_name = mutation_info.get("mutation_name")
    has_markers = any(label.endswith("m") for chain in structures.values() for label in chain)
    units: dict[str, str] = {}
    if has_markers:
        if not isinstance(mutation_name, str) or mutation_name not in rd.mutations:
            raise ValueError(f"unknown mutation {mutation_name!r}")
        unit_path = _one(root, ".mutated_unit.fa")
        paths.append(unit_path)
        units = fasta_records(unit_path)
    elif mutation_name:
        raise ValueError("mutation metadata without a marked repeat")
    warnings = []
    provenance: dict[str, Any] = dict(stats.get("provenance") or {})
    if provenance.get("seed") is None:
        warnings.append("missing_seed")
    _read_provenance(root, provenance, paths)
    provenance["vntr_coverage"] = stats.get("vntr_coverage", {})
    manifests = sorted(root.glob("*_read_truth.tsv.gz"))
    provenance["read_source_truth"] = "available" if len(manifests) == 1 else "unavailable"
    if len(manifests) == 1:
        paths.append(manifests[0])
    warnings.append("nominal_coverage_is_not_retained_molecules")
    haplotypes = []
    used_units: set[str] = set()
    targets = []
    for i, name in enumerate(expected_ids):
        labels, sequences, events = [], [], []
        details = []
        for position, label in enumerate(structures[name], 1):
            mutant = label.endswith("m")
            parent = label[:-1] if mutant else label
            if parent not in rd.repeats:
                raise ValueError(f"{name}: unknown repeat {parent!r}")
            sequence = rd.repeats[parent]
            if mutant:
                definition = rd.mutations[str(mutation_name)]
                if parent not in definition["allowed_repeats"]:
                    raise ValueError(f"{name}: mutation disallowed on {parent}")
                key = f"{name}_repeat_{position}"
                sequence = units.get(key, "")
                if sequence != _apply_mutation(rd.repeats[parent], definition["changes"]):
                    raise ValueError(f"{key}: mutated unit inconsistent with dictionary edits")
                used_units.add(key)
                events.append(Event(position, parent, mutation_name))
                details.append({"position": position, "repeat": parent})
                targets.append(f"{i + 1},{position}")
                label = f"{parent}:{mutation_name}"
            labels.append(label)
            sequences.append(sequence)
        sequence = "".join(sequences)
        if records[name] != rd.flanking_left + sequence + rd.flanking_right:
            raise ValueError(f"{name}: reconstructed sequence/flanks do not match FASTA")
        row = hap_stats[i]
        if row["repeat_count"] != len(labels) or row["vntr_length"] != len(sequence):
            raise ValueError(f"{name}: repeat count or VNTR length inconsistent")
        if row["mutant_repeat_count"] != len(events) or row["mutation_details"] != details:
            raise ValueError(f"{name}: mutation marker/statistics discrepancy")
        if row.get("snp_count", 0) or row.get("applied_snps", []):
            raise ValueError(f"{name}: SNP truth adapter not implemented")
        if row.get("repeat_lengths") != list(map(len, sequences)):
            warnings.append(f"stale_repeat_lengths:{name}")
        haplotypes.append(TruthHaplotype(name, sequence, tuple(labels), tuple(events)))
    if used_units != set(units):
        raise ValueError("unmatched mutated-unit FASTA records")
    if sorted(mutation_info.get("mutation_targets", [])) != sorted(targets):
        raise ValueError("mutation target metadata disagrees with structure")
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    return TruthSample(root.name, tuple(haplotypes), provenance, hashes, tuple(warnings))


def load_truth(sample_dir: Path, repeat_dict: RepeatDictionary) -> TruthSample:
    """Validate all authoritative artifacts before returning scorable truth."""
    try:
        return _load(sample_dir, repeat_dict)
    except (OSError, ValueError, KeyError, TypeError, IndexError, AttributeError) as exc:
        raise TruthValidationError(f"{sample_dir}: {exc}") from exc
