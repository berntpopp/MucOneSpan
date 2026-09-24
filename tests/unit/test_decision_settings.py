"""Tests for resolving clinical-decision thresholds: explicit > recorded > default."""

from __future__ import annotations

from typing import Any

import pytest

from muc_one_span import clinical_gates
from muc_one_span.decision_settings import resolve_decision_settings
from muc_one_span.settings import DEFAULT_SETTINGS, ClinicalDecisionSettings


def test_default_thresholds_and_legacy_alias() -> None:
    settings = ClinicalDecisionSettings()
    assert (settings.max_ambiguous_bases, settings.legacy_min_total_reads) == (10, 30)
    assert clinical_gates.LEGACY_MIN_TOTAL_READS == 30


def test_empty_summary_resolves_to_default() -> None:
    settings, source = resolve_decision_settings({})
    assert settings == DEFAULT_SETTINGS.clinical_decision
    assert source == "default"


def test_v0160_style_configuration_without_clinical_decision_resolves_to_default() -> None:
    """A v0.16.0-style summary carries configuration.settings but no clinical_decision."""
    summary = {"configuration": {"settings": {"run": {"threads": 8}}}}
    settings, source = resolve_decision_settings(summary)
    assert settings == DEFAULT_SETTINGS.clinical_decision
    assert source == "default"


def test_recorded_clinical_decision_is_used_when_no_explicit_override() -> None:
    summary = {
        "configuration": {
            "settings": {
                "clinical_decision": {"max_ambiguous_bases": 20, "legacy_min_total_reads": 45}
            }
        }
    }
    settings, source = resolve_decision_settings(summary)
    assert settings == ClinicalDecisionSettings(max_ambiguous_bases=20, legacy_min_total_reads=45)
    assert source == "recorded_configuration"


def test_explicit_settings_take_precedence_over_recorded() -> None:
    summary = {"configuration": {"settings": {"clinical_decision": {"max_ambiguous_bases": 20}}}}
    explicit = ClinicalDecisionSettings(max_ambiguous_bases=5)
    settings, source = resolve_decision_settings(summary, explicit)
    assert settings is explicit
    assert source == "explicit"


@pytest.mark.parametrize(
    "recorded",
    [
        {"max_ambiguous_bases": -1},
        {"unknown": 1},
    ],
)
def test_invalid_recorded_clinical_decision_raises(recorded: dict[str, Any]) -> None:
    summary = {"configuration": {"settings": {"clinical_decision": recorded}}}
    with pytest.raises(ValueError):
        resolve_decision_settings(summary)
