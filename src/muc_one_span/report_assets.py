"""Vendored assets and cryptographic verification for offline report generation."""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

#: The pinned igv.js release.
IGV_VERSION = "3.0.2"

#: Where the vendored bytes were fetched from upstream.
IGV_SOURCE_URL = "https://cdn.jsdelivr.net/npm/igv@3.0.2/dist/igv.min.js"

#: SHA-256 of the decompressed library (1,310,337 bytes of igv.min.js).
IGV_SHA256 = "ab1aa79c514ee3a0d66a0ffc788b6d37803910e62cf6d114d9b2909d96b5e790"

#: SHA-256 of the gzipped asset shipped in this package (372,690 bytes).
IGV_GZIP_SHA256 = "0d8b512654b2ef588009453c403c8f1329dce88eedca90ba9e60888af6b2f79f"

#: Directory containing vendored static assets.
ASSET_DIR = Path(__file__).resolve().parent / "assets"

#: The controlled template passed to igv-reports.
IGV_REPORT_TEMPLATE_PATH = (
    Path(__file__).resolve().parent / "templates" / "igv_report_template.html"
)

#: The library placeholder in the controlled template.
IGV_REPORT_LIBRARY_MARKER = "@MUCONESPAN_IGV_LIBRARY@"

#: The gzipped library asset path.
IGV_ASSET_PATH = ASSET_DIR / f"igv-{IGV_VERSION}.min.js.gz"

#: Markers within igv-reports generated HTML.
IGV_TABLE_JSON_MARKER = "const tableJson = "
IGV_SESSION_DICTIONARY_MARKER = "const sessionDictionary = "
IGV_CONTAINER_MARKER = '<div id="container"'
IGV_BODY_END_MARKER = "</body>"

EMPTY_TABLE_JSON = "[]"
EMPTY_SESSION_DICTIONARY = "{}"

REPORT_IGV_EMBEDDED = "embedded"
REPORT_IGV_SIDECAR = "sidecar"
REPORT_IGV_OFF = "off"
REPORT_IGV_MODES = (REPORT_IGV_EMBEDDED, REPORT_IGV_SIDECAR, REPORT_IGV_OFF)
DEFAULT_REPORT_IGV = REPORT_IGV_EMBEDDED


def _verify_bytes(payload: bytes, expected_sha256: str, description: str) -> None:
    """Validate that payload matches the expected SHA-256 digest."""
    actual = hashlib.sha256(payload).hexdigest()
    if actual != expected_sha256:
        msg = (
            f"{description} does not match pinned SHA-256: expected {expected_sha256}, "
            f"got {actual} over {len(payload)} bytes."
        )
        logger.error(msg)
        raise ValueError(msg)
    logger.debug("%s matches pinned SHA-256 (%d bytes).", description, len(payload))


def verify_asset(path: Path, expected_sha256: str) -> None:
    """Verify that a file exists and matches its expected SHA-256 digest."""
    if not path.is_file():
        msg = f"Vendored asset {path} is missing."
        logger.error(msg)
        raise ValueError(msg)
    _verify_bytes(path.read_bytes(), expected_sha256, f"Vendored asset {path.name}")


def _read_verified_asset() -> tuple[bytes, bytes]:
    """Read the vendored asset and verify both compressed and decompressed forms."""
    if not IGV_ASSET_PATH.is_file():
        msg = f"Vendored asset {IGV_ASSET_PATH} is missing."
        logger.error(msg)
        raise ValueError(msg)
    compressed = IGV_ASSET_PATH.read_bytes()
    _verify_bytes(compressed, IGV_GZIP_SHA256, f"Vendored asset {IGV_ASSET_PATH.name}")
    source = gzip.decompress(compressed)
    _verify_bytes(source, IGV_SHA256, f"Decompressed igv.js {IGV_VERSION}")
    if b"</script" in source.lower():
        msg = (
            f"Decompressed igv.js {IGV_VERSION} contains script-closing sequence and "
            "cannot be safely inserted into template."
        )
        logger.error(msg)
        raise ValueError(msg)
    return compressed, source


def igv_library_source() -> bytes:
    """Return the verified, decompressed igv.js source bytes."""
    return _read_verified_asset()[1]


def igv_payload(mode: str) -> str | None:
    """Return base64 encoded gzipped igv.js for embedded mode, or None otherwise."""
    if mode not in REPORT_IGV_MODES:
        raise ValueError(f"Unknown --report-igv mode {mode!r}; expected one of {REPORT_IGV_MODES}")
    if mode != REPORT_IGV_EMBEDDED:
        return None

    compressed, source = _read_verified_asset()
    payload = base64.b64encode(compressed).decode("ascii")
    logger.info(
        "Embedding igv.js %s: %d source bytes -> %d base64 chars (%.1f%%).",
        IGV_VERSION,
        len(source),
        len(payload),
        100.0 * len(payload) / len(source),
    )
    return payload


def igv_provenance(mode: str) -> str:
    """Return provenance string for report footer describing alignment browser status."""
    if mode not in REPORT_IGV_MODES:
        raise ValueError(f"Unknown --report-igv mode {mode!r}; expected one of {REPORT_IGV_MODES}")
    if mode == REPORT_IGV_OFF:
        return "not included (--report-igv off)"
    where = (
        "embedded in this file, gzipped"
        if mode == REPORT_IGV_EMBEDDED
        else "in separate HTML sidecar (--report-igv sidecar)"
    )
    return f"igv.js {IGV_VERSION} · sha256 {IGV_SHA256} · {where}"


def extract_line_after(content: str, marker: str) -> str:
    """Extract trimmed line text immediately following the marker."""
    start = content.find(marker)
    if start == -1:
        return ""
    start += len(marker)
    end = content.find("\n", start)
    if end == -1:
        end = len(content)
    return content[start:end].strip()


def js_json_literal(fragment: str, fallback: str) -> str:
    """Re-serialise JSON fragment escaping all '<' characters for HTML script safety."""
    candidate = fragment.strip().removesuffix(";").strip()
    if not candidate:
        return fallback
    try:
        value = json.loads(candidate)
    except ValueError as e:
        logger.warning("IGV fragment could not be parsed as JSON: %s", e)
        return fallback

    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return encoded.replace("<", "\\u003c")


def extract_igv_fragments(content: str) -> tuple[str, str, str]:
    """Extract container markup, tableJson, and sessionDictionary from generated IGV page."""
    igv_start = content.find(IGV_CONTAINER_MARKER)
    if igv_start == -1:
        logger.error("Failed to extract IGV content from report.")
        return "", "", ""

    candidates = [
        content.find("</main>", igv_start),
        content.find("<script", igv_start),
        content.find(IGV_BODY_END_MARKER, igv_start),
    ]
    valid_ends = [c for c in candidates if c != -1]
    if not valid_ends:
        logger.error("Failed to extract IGV content from report.")
        return "", "", ""

    igv_end = min(valid_ends)
    igv_content = content[igv_start:igv_end].strip()
    table_json = extract_line_after(content, IGV_TABLE_JSON_MARKER)
    session_dictionary = extract_line_after(content, IGV_SESSION_DICTIONARY_MARKER)

    return igv_content, table_json, session_dictionary
