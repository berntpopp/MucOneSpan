"""Resolve the clinical-decision thresholds used by ``report.compute_clinical_decision``.

Precedence: an explicit ``ClinicalDecisionSettings`` argument; otherwise the
section recorded at ``summary["configuration"]["settings"]["clinical_decision"]``
by a prior ``run`` (see ``cli_settings.write_run_configuration``); otherwise the
central defaults. A recorded section is validated the same way as a loaded
configuration file, so an unrecognised field or an out-of-range value raises
``ValueError`` instead of being silently ignored or defaulted.
"""

from __future__ import annotations

from typing import Any

from muc_one_span.settings import (
    DEFAULT_SETTINGS,
    ClinicalDecisionSettings,
    build_settings_section,
)


def resolve_decision_settings(
    summary: dict[str, Any], explicit: ClinicalDecisionSettings | None = None
) -> tuple[ClinicalDecisionSettings, str]:
    """Resolve the thresholds to apply and where they came from.

    Source is ``"explicit"``, ``"recorded_configuration"``, or ``"default"``.
    """
    if explicit is not None:
        return explicit, "explicit"
    configuration = summary.get("configuration")
    settings = configuration.get("settings") if isinstance(configuration, dict) else None
    recorded = settings.get("clinical_decision") if isinstance(settings, dict) else None
    if recorded is None:
        return DEFAULT_SETTINGS.clinical_decision, "default"
    section: ClinicalDecisionSettings = build_settings_section(
        "clinical_decision", ClinicalDecisionSettings, recorded
    )
    return section, "recorded_configuration"
