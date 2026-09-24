"""Calibration grids over runtime settings: expansion, strict validation, overlays.

A grid file is a JSON object mapping dotted ``section.field`` settings keys
(any `RuntimeSettings` section field, for example ``hybrid.het_af_min`` or a
ladder ``calling`` threshold) to either a non-empty list of values or an
inclusive ``{"min", "max", "step"}`` range. The Cartesian product over the
keys, in sorted key order, gives the grid points; an empty grid is the single
base point.

Each point becomes a complete schema-one settings overlay: the effective base
settings (`settings_as_dict` of the strictly loaded base file, so relative
resource paths are already absolute) with the point's values and the
calibrated engine applied. Every overlay is validated by the strict
`muc_one_span.settings.load_settings` before any run, and is addressed by the
SHA-256 of its canonical JSON (sorted keys, compact separators, as
`BenchConfig.sha256`), so equal overlays share one result directory.
"""

from __future__ import annotations

import copy
import hashlib
import itertools
import json
import tempfile
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Any

from muc_one_span.settings import load_settings, settings_as_dict

# Run options the benchmark harness passes explicitly on every case's command
# line (they override a config file's values), so a grid over them would be
# silently ignored. ``run.engine`` comes from ``calibrate --engine``.
HARNESS_OWNED = frozenset({"run.engine", "run.platform", "run.threads", "run.clair3_model"})
RANGE_KEYS = frozenset({"min", "max", "step"})
KEY_PARTS = 2  # ``section.field``
OVERLAY_FILE = "config.json"

# ``benchsim calibrate --stage``: "full" runs the whole pipeline (`run_split` +
# `evaluate`, any `RuntimeSettings` field); "lengths" fits only the hybrid length
# model (Task 15d) and accepts only the settings that model reads.
DEFAULT_STAGE = "full"
LENGTHS_STAGE = "lengths"
STAGES = (DEFAULT_STAGE, LENGTHS_STAGE)

# The exact `HybridSettings` fields read by `hybrid.spans`, `hybrid.lengths` and
# `hybrid.smear` (S1 anchor search, S2 length model, smear significance test):
# everything else cannot change a length-stage result, so a grid over it is
# refused for ``--stage lengths`` before any point is run.
LENGTH_STAGE_KEYS = frozenset(
    f"hybrid.{name}"
    for name in (
        "anchor_max_edits",
        "min_span_units",
        "max_span_units",
        "flank_anchor_bp",
        "flank_anchor_edit_divisor",
        "flank_anchor_edit_floor",
        "peak_window_base_bp",
        "peak_window_per_unit_bp",
        "kde_bandwidth_base_bp",
        "kde_bandwidth_per_bp",
        "kde_kernel_truncation_bw",
        "kde_grid_step_bp",
        "kde_grid_margin_bp",
        "peak_min_separation_units",
        "rejected_peak_noise_reads",
        "smear_short_product_units",
        "peak_far_near_boundary_units",
        "far_peak_min_frac",
        "near_peak_min_frac",
        "min_peak_reads",
        "smear_test_window_frac",
        "smear_background_min_reads",
        "smear_background_flank_units",
        "smear_test_correction",
        "smear_test_alpha",
        "smear_test_borderline_factor",
    )
)


@dataclass(frozen=True)
class GridPoint:
    """One validated grid point: its values, full overlay and content hash."""

    sha256: str
    values: dict[str, Any]
    config: dict[str, Any]


def canonical_sha256(data: Any) -> str:
    """SHA-256 of canonical JSON (sorted keys, compact separators)."""
    text = json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode()).hexdigest()


def _unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in pairs:
        if key in out:
            raise ValueError(f"duplicate grid key: {key}")
        out[key] = value
    return out


def _nonfinite(value: str) -> None:
    raise ValueError(f"Nonfinite JSON number: {value}")


def read_strict_json(path: Path) -> Any:
    """JSON with duplicate keys and nonfinite numbers rejected."""
    return json.loads(Path(path).read_text(), object_pairs_hook=_unique, parse_constant=_nonfinite)


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def expand_values(key: str, spec: Any) -> list[Any]:
    """A key's values: a list as given, or an inclusive ``min``/``max``/``step`` range.

    Ranges use exact decimal arithmetic (no accumulated float error) and stay
    integers when all three bounds are integers.
    """
    if isinstance(spec, list):
        if not spec:
            raise ValueError(f"{key}: needs at least one value")
        texts = [json.dumps(v, sort_keys=True) for v in spec]
        if len(set(texts)) != len(texts):
            raise ValueError(f"{key}: duplicate values")
        return list(spec)
    if not isinstance(spec, dict):
        raise ValueError(f"{key}: must be a list of values or a min/max/step range")
    if spec.keys() != RANGE_KEYS:
        raise ValueError(f"{key}: a range needs exactly min, max and step")
    if not all(_is_number(spec[k]) for k in RANGE_KEYS):
        raise ValueError(f"{key}: range bounds must be numbers")
    low, high, step = (Decimal(str(spec[k])) for k in ("min", "max", "step"))
    if step <= 0:
        raise ValueError(f"{key}: step must be > 0")
    if high < low:
        raise ValueError(f"{key}: max must be >= min")
    count = int((high - low) // step) + 1
    decimals = [low + i * step for i in range(count)]
    if all(isinstance(spec[k], int) for k in RANGE_KEYS):
        return [int(d) for d in decimals]
    return [float(d) for d in decimals]


def _check_key(key: str) -> None:
    parts = key.split(".")
    if len(parts) != KEY_PARTS or not all(parts):
        raise ValueError(f"grid key {key!r} must be a dotted section.field settings key")
    if key in HARNESS_OWNED:
        raise ValueError(f"grid key {key!r} is set by the harness on every run")


def load_grid(path: Path) -> dict[str, list[Any]]:
    """Read and expand a grid file (keys checked; values validated later per point)."""
    data = read_strict_json(path)
    if not isinstance(data, dict):
        raise ValueError("a calibration grid must be a JSON object")
    grid: dict[str, list[Any]] = {}
    for key in sorted(data):
        _check_key(key)
        grid[key] = expand_values(key, data[key])
    return grid


def check_stage_keys(grid: dict[str, list[Any]], stage: str) -> None:
    """Refuse a grid key that cannot affect ``stage`` before any point is run.

    The full stage allows any ``RuntimeSettings`` field (unchanged 15b/15c
    behaviour); the lengths stage allows only `LENGTH_STAGE_KEYS`.
    """
    if stage != LENGTHS_STAGE:
        return
    bad = sorted(k for k in grid if k not in LENGTH_STAGE_KEYS)
    if bad:
        raise ValueError(
            "--stage lengths calibrates only the hybrid length model "
            f"(spans/lengths/smear settings); not a length-model key: {', '.join(bad)}"
        )


def grid_points(grid: dict[str, list[Any]]) -> list[dict[str, Any]]:
    """Cartesian product of the grid's values, in sorted key order."""
    keys = sorted(grid)
    return [
        dict(zip(keys, combo, strict=True)) for combo in itertools.product(*(grid[k] for k in keys))
    ]


def base_settings(path: Path | None) -> dict[str, Any]:
    """Effective base settings (defaults, or the strictly loaded file) as a JSON dict."""
    try:
        return settings_as_dict(load_settings(Path(path) if path is not None else None))
    except TypeError as exc:
        raise ValueError(str(exc)) from exc


def overlay(base: dict[str, Any], values: dict[str, Any], engine: str) -> dict[str, Any]:
    """Base settings with ``values`` and ``run.engine`` applied (base left untouched)."""
    data = copy.deepcopy(base)
    for key, value in {**values, "run.engine": engine}.items():
        section, name = key.split(".")
        target = data.setdefault(section, {})
        if not isinstance(target, dict):
            raise ValueError(f"{key}: {section} is not a settings section")
        target[name] = value
    return data


def validate_settings(data: dict[str, Any]) -> None:
    """Validate an overlay with the strict settings loader (via a temporary file)."""
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / OVERLAY_FILE
        path.write_text(json.dumps(data, allow_nan=False))
        try:
            load_settings(path)
        except TypeError as exc:
            raise ValueError(str(exc)) from exc


def build_points(grid: dict[str, list[Any]], base: dict[str, Any], engine: str) -> list[GridPoint]:
    """Every grid point as a validated overlay; the first invalid point raises."""
    points = []
    for values in grid_points(grid):
        config = overlay(base, values, engine)
        try:
            validate_settings(config)
        except ValueError as exc:
            raise ValueError(f"grid point {json.dumps(values, sort_keys=True)}: {exc}") from exc
        points.append(GridPoint(canonical_sha256(config), values, config))
    return points
