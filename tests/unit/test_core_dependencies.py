"""The hybrid engine's libraries are core dependencies (v0.17.0); the extra stays an alias."""

from __future__ import annotations

import importlib.metadata
import re

from muc_one_span.hybrid import align, poa

HYBRID_LIBRARIES = ("edlib", "pyabpoa", "pyspoa")


def _requirements() -> list[str]:
    return importlib.metadata.requires("muc_one_span") or []


def _name(requirement: str) -> str:
    return re.split(r"[\s;<>=!~\[]", requirement, maxsplit=1)[0].lower()


def test_hybrid_libraries_are_core_dependencies() -> None:
    core = {_name(r) for r in _requirements() if "extra ==" not in r}
    assert set(HYBRID_LIBRARIES) <= core


def test_hybrid_extra_is_kept_for_backwards_compatible_installs() -> None:
    extras = importlib.metadata.metadata("muc_one_span").get_all("Provides-Extra") or []
    assert "hybrid" in extras


def test_hybrid_modules_bind_their_libraries_at_import() -> None:
    assert align.edlib.__name__ == "edlib"
    assert (poa.pyabpoa.__name__, poa.spoa.__name__) == ("pyabpoa", "spoa")
