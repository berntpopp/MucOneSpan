"""Public API defaults track central settings while preserving concrete preset defaults."""

from __future__ import annotations

import inspect
import runpy
import sys
from dataclasses import replace
from pathlib import Path
from types import ModuleType

import pytest

from muc_one_span import alleles, calling, mapping, settings

DEFAULT_FIELDS = {
    "mapping": {"map_reads": {"threads": "threads"}},
    "calling": {
        "_extract_and_remap_reads": {"threads": "threads"},
        "run_clair3": {
            "model_path": "clair3_model",
            "platform": "platform",
            "threads": "threads",
        },
        "disambiguate_same_length_alleles": {
            "clair3_model": "clair3_model",
            "threads": "threads",
            "min_qual": "min_qual",
            "platform": "platform",
        },
        "call_variants_per_allele": {
            "clair3_model": "clair3_model",
            "threads": "threads",
            "min_qual": "min_qual",
            "platform": "platform",
        },
    },
    "alleles": {"detect_alleles": {"min_coverage": "min_coverage"}},
}


def assert_default_contract(namespace: dict, module: str, run: settings.RunSettings) -> None:
    for function, fields in DEFAULT_FIELDS[module].items():
        parameters = inspect.signature(namespace[function]).parameters
        for argument, setting in fields.items():
            assert parameters[argument].default == getattr(run, setting), (function, argument)


def test_installed_api_defaults_match_central_settings() -> None:
    for module in (mapping, calling, alleles):
        assert_default_contract(
            vars(module), module.__name__.rsplit(".", 1)[1], settings.DEFAULT_SETTINGS.run
        )
    assert inspect.signature(calling.call_variants_per_allele).parameters["min_qual"].default == 5.0
    assert (
        type(inspect.signature(calling.call_variants_per_allele).parameters["min_qual"].default)
        is float
    )
    for function in (mapping.map_reads, calling.call_variants_per_allele):
        assert inspect.signature(function).parameters["preset"].default == "map-hifi"


@pytest.mark.parametrize(("preset", "expected"), [(None, "lr:hq"), ("map-pb", "map-pb")])
def test_api_imports_use_central_defaults_instead_of_duplicate_literals(
    monkeypatch: pytest.MonkeyPatch,
    preset: str | None,
    expected: str,
) -> None:
    # Execute fresh module definitions without reloading live modules or changing runtime state.
    # This simulates a future reviewed change to central defaults and catches disconnected copies.
    configured = replace(
        settings.DEFAULT_SETTINGS,
        run=replace(
            settings.DEFAULT_SETTINGS.run,
            threads=7,
            min_coverage=12,
            min_qual=8.5,
            platform="ont",
            clair3_model="reviewed_model",
            minimap2_preset=preset,
        ),
    )
    monkeypatch.setattr(settings, "DEFAULT_SETTINGS", configured)
    namespaces = {"mapping": runpy.run_path(str(Path(mapping.__file__)))}
    isolated_mapping = ModuleType("muc_one_span.mapping")
    isolated_mapping.__dict__.update(namespaces["mapping"])
    monkeypatch.setitem(sys.modules, "muc_one_span.mapping", isolated_mapping)
    namespaces["calling"] = runpy.run_path(str(Path(calling.__file__)))
    namespaces["alleles"] = runpy.run_path(str(Path(alleles.__file__)))
    for module, namespace in namespaces.items():
        assert_default_contract(namespace, module, configured.run)
    for module, functions in {
        "mapping": ("map_reads", "_run_mapping_pipeline"),
        "calling": (
            "_extract_and_remap_reads",
            "disambiguate_same_length_alleles",
            "call_variants_per_allele",
        ),
    }.items():
        for function in functions:
            assert (
                inspect.signature(namespaces[module][function]).parameters["preset"].default
                == expected
            )
