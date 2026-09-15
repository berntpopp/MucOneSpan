# src/muc_one_span/config.py
"""Configuration module for repeat dictionary loading and classification."""

from __future__ import annotations

import importlib.resources
import json
from dataclasses import dataclass
from pathlib import Path


def _bundled_repeats_path() -> Path:
    """Resolve bundled repeat dictionary using importlib.resources.

    Works both in editable installs and when installed as a package
    (e.g., from Docker or pip install).
    """
    ref = importlib.resources.files("muc_one_span.data.repeats").joinpath("repeats.json")
    # as_posix() returns a Traversable; for Path compatibility use the context
    return Path(str(ref))


@dataclass
class RepeatDictionary:
    """Holds the MUC1-VNTR repeat type definitions and flanking sequences.

    Attributes:
        repeats: Mapping of repeat ID to 60-bp consensus sequence.
        repeat_length_bp: Expected length in base pairs for each repeat unit.
        pre_repeat_ids: Ordered IDs for repeats that precede the main VNTR array.
        after_repeat_ids: Ordered IDs for repeats that follow the main VNTR array.
        canonical_repeat: ID of the canonical repeat type (usually "X").
        flanking_left: Left flanking sequence (hg38).
        flanking_right: Right flanking sequence (hg38).
        vntr_region: Genomic coordinates of the VNTR region (hg38).
        source: Provenance string for the repeat definitions.
    """

    repeats: dict[str, str]
    repeat_length_bp: int
    pre_repeat_ids: list[str]
    after_repeat_ids: list[str]
    canonical_repeat: str
    flanking_left: str
    flanking_right: str
    vntr_region: str
    source: str
    mutations: dict[str, dict]
    mutated_sequences: dict[str, tuple[str, str]]
    seq_to_id: dict[str, str]


def _apply_mutation(sequence: str, changes: list[dict]) -> str:
    """Apply mutation changes to a repeat sequence.

    Uses 1-based indexing (matching MucOneUp/Vrbacka conventions):
    ``start=60`` means "insert before position 60 (1-based)" which is
    0-based index 59.  Python's ``list.insert(i, x)`` inserts before
    index ``i``, so we use ``insert(start - 1, ...)``.

    Changes are applied in reverse order to preserve indices.
    """
    result = list(sequence)
    for change in sorted(changes, key=lambda c: c["start"], reverse=True):
        start = change["start"] - 1  # convert to 0-based
        if change["type"] == "insert":
            result[start:start] = list(change["sequence"])
        elif change["type"] == "delete":
            end = change["end"]  # 1-based inclusive end == 0-based exclusive end
            del result[start:end]
        elif change["type"] == "delete_insert":
            # In MucOneUp/Vrbacka conventions, start and end are 1-based retained boundaries.
            # Bases strictly between start and end are deleted, and sequence is inserted.
            end_idx = change["end"] - 1
            del result[start + 1 : end_idx]
            result[start + 1 : start + 1] = list(change["sequence"])
    return "".join(result)


def _precompute_mutated_sequences(
    repeats: dict[str, str],
    mutations: dict[str, dict],
) -> dict[str, tuple[str, str]]:
    """Pre-compute mutated sequences for all (repeat, mutation) combinations.

    Returns:
        Dict mapping mutated DNA sequence -> (parent_repeat_id, mutation_name).
    """
    result: dict[str, tuple[str, str]] = {}
    for mut_name, mut_def in mutations.items():
        for repeat_id in mut_def["allowed_repeats"]:
            if repeat_id not in repeats:
                continue
            seq = repeats[repeat_id]
            mutated = _apply_mutation(seq, mut_def["changes"])
            result[mutated] = (repeat_id, mut_name)
    return result


def load_repeat_dictionary(path: Path | None = None) -> RepeatDictionary:
    """Load the repeat dictionary from a JSON file.

    Args:
        path: Path to a JSON file.  If *None*, the bundled
              ``data/repeats/repeats.json`` is used.

    Returns:
        A :class:`RepeatDictionary` instance populated from the file.

    Raises:
        FileNotFoundError: If the JSON file does not exist.
        KeyError: If required keys are missing from the JSON.
    """
    resolved = path if path is not None else _bundled_repeats_path()
    if not resolved.exists():
        raise FileNotFoundError(f"Repeat dictionary not found: {resolved}")

    data = json.loads(resolved.read_text(encoding="utf-8"))

    mutations_raw = data.get("mutations", {})
    # Filter to real mutations only (exclude any with type=synthetic)
    mutations = {k: v for k, v in mutations_raw.items() if v.get("type", "real") != "synthetic"}
    mutated_seqs = _precompute_mutated_sequences(data["repeats"], mutations)

    flanking = data["flanking_hg38"]
    return RepeatDictionary(
        repeats=data["repeats"],
        repeat_length_bp=data["repeat_length_bp"],
        pre_repeat_ids=data["pre_repeat_ids"],
        after_repeat_ids=data["after_repeat_ids"],
        canonical_repeat=data["canonical_repeat"],
        flanking_left=flanking["left"],
        flanking_right=flanking["right"],
        vntr_region=flanking["vntr_region"],
        source=data["source"],
        mutations=mutations,
        mutated_sequences=mutated_seqs,
        seq_to_id={seq: rid for rid, seq in data["repeats"].items()},
    )


def classify_repeat_id(repeat_id: str, dictionary: RepeatDictionary | None = None) -> str:
    """Classify a repeat ID into a structural category.

    Categories:
    - ``"pre"``       - precedes the main VNTR array (e.g. "1" to "5C")
    - ``"after"``     - follows the main VNTR array (e.g. "6" to "9")
    - ``"canonical"`` - the canonical repeat type (usually "X")
    - ``"variable"``  - any other known repeat variant (e.g. "A", "B", ...)

    Args:
        repeat_id: The repeat identifier to classify.
        dictionary: Optional :class:`RepeatDictionary` to use.  If *None*, the
                    bundled dictionary is loaded automatically.

    Returns:
        One of ``"pre"``, ``"after"``, ``"canonical"``, or ``"variable"``.
    """
    if dictionary is None:
        dictionary = load_repeat_dictionary()

    if repeat_id in dictionary.pre_repeat_ids:
        return "pre"
    if repeat_id in dictionary.after_repeat_ids:
        return "after"
    if repeat_id == dictionary.canonical_repeat:
        return "canonical"
    return "variable"
