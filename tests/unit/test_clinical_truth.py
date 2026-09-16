"""Independent sequence-truth acquisition and checksum contracts."""

import hashlib
import io
import json
from pathlib import Path

import pytest

from muc_one_span.clinical_truth import hydrate_truth


def test_sequence_truth_is_reverse_complemented_once_and_checked(tmp_path: Path, monkeypatch):
    sequence = "AGTC"
    source = {
        "url": "https://example/sequence",
        "strand": "-",
        "sha256": hashlib.sha256(b"GACT").hexdigest(),
        "length": 4,
    }
    ledger = {"samples": [{"sequence_truth": {"sequence_records": [source]}}]}
    monkeypatch.setattr(
        "muc_one_span.clinical_truth.urllib.request.urlopen",
        lambda *a, **k: io.BytesIO(json.dumps({"dna": sequence}).encode()),
    )
    result = hydrate_truth(ledger, tmp_path)
    assert result["samples"][0]["sequence_truth"]["sequences"] == ["GACT"]
    monkeypatch.setattr(
        "muc_one_span.clinical_truth.urllib.request.urlopen",
        lambda *a, **k: pytest.fail("must replay validated cache offline"),
    )
    assert hydrate_truth(ledger, tmp_path) == result
    assert "sequences" not in ledger["samples"][0]["sequence_truth"]


def test_sequence_checksum_mismatch_is_not_published(tmp_path: Path, monkeypatch):
    ledger = {
        "samples": [
            {
                "sequence_truth": {
                    "sequence_records": [
                        {"url": "https://example/s", "strand": "+", "sha256": "0" * 64, "length": 4}
                    ]
                }
            }
        ]
    }
    monkeypatch.setattr(
        "muc_one_span.clinical_truth.urllib.request.urlopen",
        lambda *a, **k: io.BytesIO(b'{"dna":"AGTC"}'),
    )
    with pytest.raises(ValueError, match="checksum"):
        hydrate_truth(ledger, tmp_path)
    assert not list(tmp_path.iterdir())
