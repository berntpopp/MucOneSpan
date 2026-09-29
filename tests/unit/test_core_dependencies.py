"""The default hybrid engine's libraries are core (v0.17.0); pyspoa stays in the `hybrid` extra."""

from __future__ import annotations

import importlib.metadata
import re

from muc_one_span.hybrid import align, poa

CORE_LIBRARIES = ("edlib", "pyabpoa")
EXTRA_ONLY = "pyspoa"


def _requirements() -> list[str]:
    return importlib.metadata.requires("muc_one_span") or []


def _name(requirement: str) -> str:
    return re.split(r"[\s;<>=!~\[]", requirement, maxsplit=1)[0].lower()


def test_default_backend_libraries_are_core_dependencies() -> None:
    core = {_name(r) for r in _requirements() if "extra ==" not in r}
    assert set(CORE_LIBRARIES) <= core
    assert EXTRA_ONLY not in core


def test_pyspoa_is_only_in_the_hybrid_extra() -> None:
    extra = {_name(r) for r in _requirements() if 'extra == "hybrid"' in r.replace("'", '"')}
    assert EXTRA_ONLY in extra


def test_hybrid_extra_is_kept_for_backwards_compatible_installs() -> None:
    extras = importlib.metadata.metadata("muc_one_span").get_all("Provides-Extra") or []
    assert "hybrid" in extras


def test_core_libraries_are_bound_at_import_and_spoa_is_lazy() -> None:
    assert align.edlib.__name__ == "edlib"
    assert poa.pyabpoa.__name__ == "pyabpoa"
    assert not hasattr(poa, "spoa")
