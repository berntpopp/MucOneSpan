"""Read caller artifacts and explicit inventories without invoking caller tools."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .models import Event, PredictedAllele, RunObservation
from .truth import fasta_records

INPUT_SUFFIXES = (".bam", ".fastq", ".fq", ".fastq.gz", ".fq.gz")


def discover_input(sample_dir: Path, selected: str | None = None) -> Path:
    """Select an explicit supported input or require exactly one candidate."""
    if selected is not None:
        path = sample_dir / selected
        if not path.is_file() or not path.name.endswith(INPUT_SUFFIXES):
            raise ValueError(f"invalid selected input: {path}")
        return path
    paths = sorted(
        p for p in sample_dir.iterdir() if p.is_file() and p.name.endswith(INPUT_SUFFIXES)
    )
    if len(paths) != 1:
        raise ValueError(f"missing or ambiguous input in {sample_dir}: {len(paths)} candidates")
    return paths[0]


def read_inventory(path: Path) -> list[dict[str, Any]]:
    """Read a nonempty JSON array of IDs or objects with sample/path metadata."""
    raw = json.loads(path.read_text())
    if not isinstance(raw, list) or not raw:
        raise ValueError("expected a nonempty sample inventory JSON array")
    rows: list[dict[str, Any]] = []
    seen = set()
    for item in raw:
        row = {"sample": item} if isinstance(item, str) else item
        if not isinstance(row, dict) or not isinstance(row.get("sample"), str):
            raise ValueError("inventory entries require a string sample identifier")
        name = row["sample"]
        if not name or name in (".", "..") or "/" in name or "\\" in name:
            raise ValueError(f"invalid sample identifier: {name!r}")
        if name in seen:
            raise ValueError(f"duplicate sample identifier: {name}")
        seen.add(name)
        rows.append(row)
    return rows


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    return value


def _prediction(
    root: Path, name: str, classification: dict[str, Any], allele: dict[str, Any]
) -> PredictedAllele:
    if "/" in name or "\\" in name or name in (".", ".."):
        raise ValueError("invalid allele identifier")
    records = fasta_records(root / f"consensus_{name}.fa")
    if len(records) != 1:
        raise ValueError(f"{name}: expected exactly one consensus record")
    structure = classification["structure"]
    mutations = classification.get("mutations", [])
    if not isinstance(structure, str) or not isinstance(mutations, list):
        raise ValueError(f"{name}: malformed classification")
    for field in (
        "sequence_source",
        "phase_status",
        "genotype_status",
        "sequence_identity_status",
        "reconstruction_status",
    ):
        value = allele.get(field, classification.get(field))
        if value is not None and not isinstance(value, str):
            raise ValueError(f"{name}: {field} must be a string")
    events = []
    for m in mutations:
        index = m.get("repeat_index")
        parent, mutation_name = m.get("closest_type"), m.get("mutation_name")
        if index is not None:
            index = _integer(index, "repeat_index")
        if any(x is not None and not isinstance(x, str) for x in (parent, mutation_name)):
            raise ValueError(f"{name}: malformed event annotation")
        if not isinstance(m.get("vcf_support_status", "unknown"), str):
            raise ValueError(f"{name}: malformed variant support status")
        read_support = m.get("read_support", {"status": "unknown"})
        if not isinstance(read_support, dict) or not isinstance(read_support.get("status"), str):
            raise ValueError(f"{name}: malformed read support")
        events.append(
            Event(
                index,
                parent,
                mutation_name,
                m.get("frameshift") is True,
                m.get("template_match") is True,
                m.get("vcf_support") is True,
                m.get("vcf_support_status", "unknown"),
                read_support["status"],
            )
        )
    from muc_one_span.settings import DEFAULT_SETTINGS

    fixed_count = _integer(
        allele.get("fixed_repeat_count", DEFAULT_SETTINGS.reference_layout.fixed_repeat_count),
        "fixed_repeat_count",
    )
    if fixed_count < 0:
        raise ValueError("fixed_repeat_count must be nonnegative")
    return PredictedAllele(
        name,
        next(iter(records.values())),
        tuple(structure.split()),
        _integer(allele["length"], "length"),
        _integer(allele["canonical_repeats"], "canonical_repeats"),
        tuple(events),
        classification.get("reconstruction_status", "unknown"),
        allele.get("sequence_source", classification.get("sequence_source")),
        allele.get(
            "independent_haplotype_evidence", classification.get("independent_haplotype_evidence")
        )
        is True,
        allele.get("phase_status", classification.get("phase_status", "unknown")),
        allele.get("genotype_status", classification.get("genotype_status", "unknown")),
        allele.get(
            "sequence_identity_status", classification.get("sequence_identity_status", "unknown")
        ),
        fixed_count,
    )


def _statuses_allow_completion(allele: dict[str, Any], classification: dict[str, Any]) -> bool:
    """Recognize producer statuses without promoting new or unknown evidence.

    Legacy artifacts may omit these additive fields. Explicit unknown/null values
    are different from absence. Check both sources so one layer cannot conceal
    another layer's unresolved state. Completed segmentation still does not imply
    empirical reference-base confidence.
    """
    accepted = {
        "phase_status": {
            "phased",
            "single_heterozygous_unordered",
            "no_informative_heterozygosity",
        },
        "reconstruction_status": {
            "complete_segmentation",
            "candidate_reference_confidence_unverified",
        },
    }
    return all(
        field not in source or (isinstance(source[field], str) and source[field] in allowed)
        for source in (allele, classification)
        for field, allowed in accepted.items()
    )


def load_observation(result_dir: Path, run_record: dict[str, Any] | None = None) -> RunObservation:
    """Preserve observations and known failures, even with a damaged status file.

    A valid typed insufficient-evidence sidecar distinguishes the CLI's expected
    nonzero coverage outcome from tool failure. An invalid sidecar cannot replace
    external evidence that execution failed or timed out.
    """
    record = dict(run_record or {})
    status = record.get("status")
    execution_failed = record.get("exit_code") not in (None, 0) or status in (
        "execution_failed",
        "timeout",
    )
    sidecar_path = result_dir / "run_status.json"
    if sidecar_path.exists():
        try:
            sidecar = json.loads(sidecar_path.read_text())
            if (
                not isinstance(sidecar, dict)
                or sidecar.get("schema_version") != 1
                or sidecar.get("status")
                not in (
                    "running",
                    "completed",
                    "execution_failed",
                    "insufficient_evidence",
                    "interrupted",
                )
            ):
                raise ValueError("invalid run_status.json schema/status")
            record["run_status"] = sidecar
            # Coverage no-calls use exit 1. A timeout or different nonzero code
            # cannot be explained by a retained insufficient-evidence sidecar.
            if status == "timeout" or record.get("exit_code") not in (None, 0, 1):
                return RunObservation(
                    "execution_failed", run_record=record, error=record.get("error")
                )
            if sidecar["status"] in ("execution_failed", "insufficient_evidence", "interrupted"):
                return RunObservation(
                    sidecar["status"], run_record=record, error=sidecar.get("error")
                )
        except (OSError, ValueError, TypeError) as exc:
            if execution_failed:
                return RunObservation(
                    "execution_failed",
                    run_record=record,
                    warnings=(f"invalid_run_status_sidecar:{exc}",),
                    error=record.get("error"),
                )
            return RunObservation("invalid_artifacts", run_record=record, error=str(exc))
    if execution_failed:
        return RunObservation("execution_failed", run_record=record, error=record.get("error"))
    if record.get("run_status", {}).get("status") == "running":
        return RunObservation(
            "invalid_artifacts", run_record=record, error="run incomplete: status is running"
        )
    if status == "not_attempted":
        return RunObservation("not_attempted", run_record=record)
    summary_path = result_dir / "summary.json"
    if not summary_path.exists():
        return RunObservation(
            "invalid_artifacts" if record.get("exit_code") == 0 else "not_attempted",
            run_record=record,
            error="missing summary.json",
        )
    warnings: tuple[str, ...] = (
        () if record.get("exit_code") == 0 else ("unknown_execution_provenance",)
    )
    try:
        data = json.loads(summary_path.read_text())
        classifications, alleles = data["classifications"], data["alleles"]
        if not isinstance(classifications, dict) or not isinstance(alleles, dict):
            raise ValueError("classifications and alleles must be objects")
        expected = {k for k, v in alleles.items() if isinstance(v, dict)}
        observed = set(classifications)
        if not observed <= expected:
            raise ValueError("classification/allele identifiers disagree")
        unresolved = expected - observed
        for name in sorted(unresolved):
            alias = alleles[name]
            target = alias.get("candidate_duplicate_of")
            if (
                alias.get("reconstruction_status") != "not_separately_resolved"
                or not isinstance(target, str)
                or target not in observed
            ):
                raise ValueError("classification/allele identifiers disagree")
            warnings += (f"unresolved_allele_alias:{name}:{target}",)
        repeats_path = result_dir / "repeats.json"
        if repeats_path.exists():
            complete = json.loads(repeats_path.read_text())
            if not isinstance(complete, dict) or set(complete) != observed:
                raise ValueError("repeats/classification identifiers disagree")
        predictions = tuple(
            _prediction(result_dir, k, v, alleles[k]) for k, v in sorted(classifications.items())
        )
        status = "completed" if predictions else "insufficient_evidence"
        if unresolved or any(
            not _statuses_allow_completion(alleles[p.name], classifications[p.name])
            or any(base not in "ACGT" for base in p.sequence)
            for p in predictions
        ):
            status = "ambiguous_reconstruction"
        return RunObservation(status, predictions, record, warnings)
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        return RunObservation(
            "invalid_artifacts", run_record=record, warnings=warnings, error=str(exc)
        )
