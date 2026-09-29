"""``muconespan settings``: print and validate schema-one runtime settings."""

from __future__ import annotations

import json
from dataclasses import fields, is_dataclass
from pathlib import Path

import click

from muc_one_span.cli_settings import current_settings
from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings, load_settings, settings_as_dict

# Sections are the dataclass-valued fields of RuntimeSettings, so a new section
# is selectable without editing this module.
SECTION_NAMES = tuple(
    f.name for f in fields(RuntimeSettings) if is_dataclass(getattr(DEFAULT_SETTINGS, f.name))
)
# Formatting of the printed JSON; it matches examples/runtime-settings.json.
JSON_INDENT = 2


def _load(path: Path) -> RuntimeSettings:
    """Run the strict loader, turning its first error into a CLI error."""
    try:
        return load_settings(path)
    except (OSError, ValueError, TypeError) as exc:
        raise click.ClickException(f"Invalid configuration {path}: {exc}") from exc


def render_settings(settings: RuntimeSettings, section: str | None = None) -> str:
    """Return settings as schema-one JSON text that ``--config`` loads back.

    With ``section``, only that section and ``schema_version`` are included; the
    omitted sections keep their central defaults when the text is loaded.
    """
    data = settings_as_dict(settings)
    if section is not None:
        data = {"schema_version": data["schema_version"], section: data[section]}
    return json.dumps(data, indent=JSON_INDENT) + "\n"


@click.group("settings")
def settings_group() -> None:
    """Print or validate runtime settings JSON."""


@settings_group.command("show")
@click.option(
    "--section",
    type=click.Choice(SECTION_NAMES),
    default=None,
    help="Print only this settings section.",
)
@click.option(
    "--config",
    "configuration",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help="Settings JSON to resolve; defaults to the global --config, else the defaults.",
)
def show(section: str | None, configuration: Path | None) -> None:
    """Print the effective settings as schema-one JSON loadable with --config."""
    settings = _load(configuration) if configuration is not None else current_settings()
    click.echo(render_settings(settings, section), nl=False)


@settings_group.command("validate")
@click.argument("file", type=click.Path(path_type=Path, dir_okay=False))
def validate(file: Path) -> None:
    """Run the strict loader on FILE; exit non-zero with the first error."""
    _load(file)
    click.echo(f"{file}: valid schema-one settings")
