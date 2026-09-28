"""Version values outside Python that must follow ``muc_one_span.version`` (ledger L258)."""

from __future__ import annotations

import re
from pathlib import Path

from muc_one_span.version import __version__

ROOT = Path(__file__).resolve().parents[2]


def test_apptainer_definition_uses_the_package_version() -> None:
    text = (ROOT / "docker" / "muconespan.def").read_text()
    image = re.search(r"^From: ghcr\.io/berntpopp/muconespan:(\S+)$", text, re.M)
    label = re.search(r"^\s+Version (\S+)$", text, re.M)
    assert image is not None and image[1] == __version__
    assert label is not None and label[1] == __version__


def test_citation_page_uses_the_package_version() -> None:
    text = (ROOT / "docs" / "about" / "citation.md").read_text()
    assert f"version = {{{__version__}}}," in text
