"""Deprecated and ignored run options: the warnings a run prints and the records it keeps.

``summary.json["deprecations"]`` is an additive list; it is empty unless a run selected a
deprecated option. v0.17.0 deprecates the ladder engine (``run.engine = "ladder"``); the
hybrid engine is the default. The ladder stays selectable until a later release removes it.

``summary.json["ignored_options"]`` (and ``run_configuration.json``) lists ladder-only
options that a hybrid run was given explicitly (command line) or through a non-default
configuration value; the hybrid engine does not use them.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import click
from click.core import ParameterSource

from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings

LADDER_ENGINE = "ladder"
LADDER_DEPRECATION_MESSAGE = (
    "The ladder engine (--engine ladder / run.engine = ladder) is deprecated and will be "
    "removed in a future release; the hybrid engine is the default. See the migration "
    "guide in the documentation."
)
# run.<field> -> CLI flag, for options only the ladder engine uses.
LADDER_ONLY_OPTIONS = {
    "clair3_model": "--clair3-model",
    "min_qual": "--min-qual",
    "minimap2_preset": "--minimap2-preset",
}
IGNORED_OPTION_WARNING = (
    "{option} is ignored by the hybrid engine; use --engine ladder (deprecated)"
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


def explicit_command_line_options() -> set[str]:
    """Ladder-only ``run`` parameters given on the command line of the current invocation."""
    ctx = click.get_current_context(silent=True)
    if ctx is None:
        return set()
    return {
        name
        for name in LADDER_ONLY_OPTIONS
        if ctx.get_parameter_source(name) == ParameterSource.COMMANDLINE
    }


def ignored_options(settings: RuntimeSettings, explicit: Iterable[str]) -> list[dict[str, Any]]:
    """Ladder-only options a non-ladder run was given explicitly or by a non-default value."""
    if settings.run.engine == LADDER_ENGINE:
        return []
    chosen = set(explicit)
    records = []
    for name, option in LADDER_ONLY_OPTIONS.items():
        value = getattr(settings.run, name)
        if name in chosen or value != getattr(DEFAULT_SETTINGS.run, name):
            records.append(
                {
                    "option": option,
                    "setting": f"run.{name}",
                    "value": value,
                    "status": "ignored",
                    "engine": settings.run.engine,
                }
            )
    return records


def warn_ignored(records: Iterable[dict[str, Any]]) -> None:
    """Print one warning line on stderr for each ignored option."""
    for record in records:
        click.echo("Warning: " + IGNORED_OPTION_WARNING.format(option=record["option"]), err=True)
