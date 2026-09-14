"""Click configuration precedence and explicit effective-run provenance."""

from __future__ import annotations

import hashlib
import json
from dataclasses import fields, replace
from pathlib import Path
from typing import Any

import click

from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings, load_settings, settings_as_dict


def configure_context(ctx: click.Context, path: Path | None) -> None:
    """Apply file defaults only; Click's explicit command-line values take precedence."""
    try:
        settings = load_settings(path)
    except (OSError, ValueError, TypeError) as exc:
        raise click.BadParameter(f"Invalid configuration: {exc}", param_hint="--config") from exc
    ctx.obj = settings
    ctx.meta["configuration_path"] = path
    values = settings_as_dict(settings)["run"]
    values["repeats_db"] = settings.repeat_dictionary
    values["flank_length"] = settings.consensus.flank_length
    values["min_units"] = settings.reference_layout.min_units
    values["max_units"] = settings.reference_layout.max_units
    if not isinstance(ctx.command, click.Group):
        raise TypeError("Configuration defaults require a Click command group")
    defaults = dict(ctx.default_map or {})
    for name, command in ctx.command.commands.items():
        configured = dict(defaults.get(name, {}))
        for parameter in command.params:
            if parameter.name in values:
                value = values[parameter.name]
                if value is not None:
                    if path is not None:
                        configured[parameter.name] = value
                    else:
                        configured.setdefault(parameter.name, value)
        defaults[name] = configured
    ctx.default_map = defaults


def current_settings() -> RuntimeSettings:
    """Return the immutable invocation settings, including for direct library use."""
    ctx = click.get_current_context(silent=True)
    return ctx.obj if ctx is not None and isinstance(ctx.obj, RuntimeSettings) else DEFAULT_SETTINGS


def current_configuration_path() -> Path | None:
    """Return the input configuration path without process-global state."""
    ctx = click.get_current_context(silent=True)
    return ctx.meta.get("configuration_path") if ctx is not None else None


def digest(path: Path) -> str:
    """Hash file bytes without reading a complete sequencing file into memory."""
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def effective_run_settings(settings: RuntimeSettings, **values: Any) -> RuntimeSettings:
    """Validate and freeze the values actually selected by Click after overrides."""
    for name in ("reference", "clair3_model"):
        if values.get(name):
            values[name] = str(Path(values[name]).resolve())
    try:
        return replace(settings, run=replace(settings.run, **values))
    except (ValueError, TypeError) as exc:
        raise click.BadParameter(str(exc)) from exc


def validate_stage_options() -> None:
    """Validate explicit stage overrides against the same contracts as JSON values."""
    ctx = click.get_current_context()
    settings = current_settings()
    values = ctx.params
    try:
        replace(
            settings.run,
            **{f.name: values[f.name] for f in fields(settings.run) if f.name in values},
        )
        replace(
            settings.reference_layout,
            **{name: values[name] for name in ("min_units", "max_units") if name in values},
        )
        if "flank_length" in values:
            replace(settings.consensus, flank_length=values["flank_length"])
    except (ValueError, TypeError) as exc:
        raise click.BadParameter(str(exc)) from exc


def write_run_configuration(
    settings: RuntimeSettings,
    configuration: Path | None,
    input_path: Path,
    reference: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Persist effective settings and input provenance before launching tools."""
    from muc_one_span.config import _bundled_repeats_path
    from muc_one_span.mapping import PLATFORM_PRESETS

    dictionary = (
        Path(settings.repeat_dictionary) if settings.repeat_dictionary else _bundled_repeats_path()
    )
    record = {
        "schema_version": 1,
        "settings": settings_as_dict(settings),
        "resolved_minimap2_preset": settings.run.minimap2_preset
        or PLATFORM_PRESETS[settings.run.platform],
        "configuration_path": str(configuration.resolve()) if configuration else None,
        "configuration_sha256": digest(configuration) if configuration else None,
        "input_path": str(input_path.resolve()),
        "input_sha256": digest(input_path),
        "resolved_reference": str(reference.resolve()),
        "reference_sha256": digest(reference),
        "resolved_repeat_dictionary": str(dictionary.resolve()),
        "repeat_dictionary_sha256": digest(dictionary),
        "precedence": "explicit CLI > configuration file > central defaults",
        "tool_selection": "executable names resolved by existing external-tool PATH handling",
        "model_selection": "explicit path"
        if settings.run.clair3_model
        else "external default unverified",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run_configuration.json").write_text(json.dumps(record, indent=2) + "\n")
    return record
