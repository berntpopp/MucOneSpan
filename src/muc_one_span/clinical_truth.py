"""Retrieve explicitly versioned sequence truth without committing participant sequences."""

from __future__ import annotations

import hashlib
import json
import urllib.request
from copy import deepcopy
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from muc_one_span.clinical_provenance import object_hash, write_json


def hydrate_truth(ledger: dict[str, Any], cache_root: Path) -> dict[str, Any]:
    """Resolve UCSC sequence descriptors into checksum-verified local truth.

    Source intervals must already be curated as zero-based half-open coordinates.
    Reverse-complement only when explicitly declared, then verify the sequence
    length and SHA256 before publishing a reusable cache entry.
    """
    result = deepcopy(ledger)
    cache_root.mkdir(parents=True, exist_ok=True)
    for sample in result["samples"]:
        truth = sample.get("sequence_truth")
        if truth is None or "sequence_records" not in truth:
            continue
        sequences = []
        for descriptor in truth["sequence_records"]:
            url = descriptor["url"]
            if urlparse(url).scheme != "https" or descriptor["strand"] not in ("+", "-"):
                raise ValueError("sequence source requires HTTPS and explicit strand")
            cache = cache_root / (object_hash(descriptor) + ".json")
            if cache.exists():
                raw = json.loads(cache.read_text())
            else:
                with urllib.request.urlopen(url, timeout=60) as response:
                    raw = json.loads(response.read(1024 * 1024))
            dna = raw["dna"].upper()
            if not dna or set(dna) - set("ACGT"):
                raise ValueError("invalid independent sequence alphabet")
            if descriptor["strand"] == "-":
                dna = dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            digest = hashlib.sha256(dna.encode("ascii")).hexdigest()
            if digest != descriptor["sha256"] or len(dna) != descriptor["length"]:
                raise ValueError("independent sequence checksum/length mismatch")
            if not cache.exists():
                write_json(cache, raw)
            sequences.append(dna)
        truth["sequences"] = sequences
        truth["sha256"] = [d["sha256"] for d in truth["sequence_records"]]
    return result
