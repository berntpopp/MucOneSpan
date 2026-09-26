"""Deprecated run options: the warning a run prints and the record ``summary.json`` keeps.

``summary.json["deprecations"]`` is an additive list; it is empty unless a run selected a
deprecated option. v0.17.0 deprecates the ladder engine (``run.engine = "ladder"``); the
hybrid engine is the default. The ladder stays selectable until a later release removes it.
"""

from __future__ import annotations

import click

from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings

LADDER_ENGINE = "ladder"
LADDER_DEPRECATION_MESSAGE = (
    "The ladder engine (--engine ladder / run.engine = ladder) is deprecated and will be "
    "removed in a future release; the hybrid engine is the default. See the migration "
    "guide in the documentation."
)


def run_deprecations(settings: RuntimeSettings) -> list[dict[str, str]]:
    """Deprecated options selected by ``settings``, one record per option."""
    if settings.run.engine != LADDER_ENGINE:
        return []
    return [
        {
            "setting": "run.engine",
            "value": LADDER_ENGINE,
            "status": "deprecated",
            "replacement": DEFAULT_SETTINGS.run.engine,
            "message": LADDER_DEPRECATION_MESSAGE,
        }
    ]


def warn_deprecations(settings: RuntimeSettings) -> None:
    """Print one warning line on stderr for each deprecated option ``settings`` selects."""
    for record in run_deprecations(settings):
        click.echo(f"Warning: {record['message']}", err=True)
