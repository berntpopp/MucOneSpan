"""Calibration grids: expansion, strict validation and content-addressed overlays."""

import ast
import inspect
import json
from pathlib import Path
from types import ModuleType

import pytest

from muc_one_span.benchsim.calibration_grid import (
    DEFAULT_STAGE,
    HARNESS_OWNED,
    LENGTH_STAGE_KEYS,
    LENGTHS_STAGE,
    STAGES,
    base_settings,
    build_points,
    check_stage_keys,
    expand_values,
    grid_points,
    load_grid,
)
from muc_one_span.hybrid import lengths as hybrid_lengths
from muc_one_span.hybrid import smear as hybrid_smear
from muc_one_span.hybrid import spans as hybrid_spans
from muc_one_span.settings import DEFAULT_SETTINGS, load_settings, settings_as_dict

# Parameter names the length-model modules use for their ``HybridSettings`` argument.
SETTINGS_PARAM_NAMES = ("settings", "h")

WINDOW = DEFAULT_SETTINGS.hybrid.smear_test_window_frac
HET = DEFAULT_SETTINGS.hybrid.het_af_min
SEED = DEFAULT_SETTINGS.hybrid.seed


def _write(tmp_path: Path, data: object, name: str = "grid.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(data))
    return path


def test_list_values_expand_to_the_cartesian_product_in_key_order() -> None:
    points = grid_points({"hybrid.seed": [SEED, SEED + 1], "hybrid.het_af_min": [HET]})
    assert points == [
        {"hybrid.het_af_min": HET, "hybrid.seed": SEED},
        {"hybrid.het_af_min": HET, "hybrid.seed": SEED + 1},
    ]


def test_empty_grid_is_the_single_base_point() -> None:
    assert grid_points({}) == [{}]


def test_integer_range_is_inclusive_and_stays_integer() -> None:
    assert expand_values("hybrid.seed", {"min": SEED, "max": SEED + 2, "step": 1}) == [
        SEED,
        SEED + 1,
        SEED + 2,
    ]


def test_float_range_has_no_accumulated_rounding() -> None:
    values = expand_values("hybrid.smear_test_window_frac", {"min": 0.15, "max": 0.35, "step": 0.1})
    assert values == [0.15, 0.25, 0.35]
    assert all(isinstance(v, float) for v in values)


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ([], "at least one value"),
        ([WINDOW, WINDOW], "duplicate"),
        ({"min": 0.1, "max": 0.3}, "exactly"),
        ({"min": 0.1, "max": 0.3, "step": 0}, "step must be > 0"),
        ({"min": 0.3, "max": 0.1, "step": 0.1}, "max must be >= min"),
        ({"min": True, "max": 1, "step": 1}, "must be numbers"),
        ("0.25", "list of values or"),
    ],
)
def test_invalid_value_specs_are_rejected(spec: object, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        expand_values("hybrid.smear_test_window_frac", spec)


def test_load_grid_rejects_non_objects_duplicates_and_bad_keys(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="JSON object"):
        load_grid(_write(tmp_path, [1]))
    dup = tmp_path / "dup.json"
    dup.write_text('{"hybrid.seed": [1], "hybrid.seed": [2]}')
    with pytest.raises(ValueError, match="duplicate"):
        load_grid(dup)
    with pytest.raises(ValueError, match=r"section\.field"):
        load_grid(_write(tmp_path, {"repeat_dictionary": ["x"]}))
    nan = tmp_path / "nan.json"
    nan.write_text('{"hybrid.seed": [NaN]}')
    with pytest.raises(ValueError, match="Nonfinite"):
        load_grid(nan)


@pytest.mark.parametrize("key", sorted(HARNESS_OWNED))
def test_harness_owned_keys_are_refused(tmp_path: Path, key: str) -> None:
    with pytest.raises(ValueError, match="set by the harness"):
        load_grid(_write(tmp_path, {key: ["x"]}))


def test_every_point_is_validated_by_the_strict_settings_loader() -> None:
    base = settings_as_dict(DEFAULT_SETTINGS)
    with pytest.raises(ValueError, match="Unknown hybrid fields"):
        build_points({"hybrid.no_such_field": [1]}, base, "hybrid")
    with pytest.raises(ValueError, match="Unknown configuration fields"):
        build_points({"nosection.field": [1]}, base, "hybrid")
    with pytest.raises(ValueError, match=r"hybrid\.smear_test_window_frac"):
        build_points({"hybrid.smear_test_window_frac": [WINDOW, 0]}, base, "hybrid")
    with pytest.raises(ValueError, match=r"run\.engine"):
        build_points({}, base, "no-such-engine")


def test_points_are_content_addressed_and_carry_the_engine() -> None:
    base = settings_as_dict(DEFAULT_SETTINGS)
    grid = {"hybrid.smear_test_window_frac": [WINDOW, WINDOW / 2]}
    first, second = build_points(grid, base, "hybrid")
    again = build_points(grid, base, "hybrid")
    assert [p.sha256 for p in again] == [first.sha256, second.sha256]
    assert first.sha256 != second.sha256
    assert first.values == {"hybrid.smear_test_window_frac": WINDOW}
    assert first.config["run"]["engine"] == "hybrid"
    assert second.config["hybrid"]["smear_test_window_frac"] == WINDOW / 2
    ladder = build_points(grid, base, "ladder")[0]
    assert ladder.sha256 != first.sha256  # the engine is part of the overlay


def test_overlay_is_loadable_and_base_paths_stay_absolute(tmp_path: Path) -> None:
    folder = tmp_path / "configs"
    folder.mkdir()
    base_file = _write(
        folder, {"schema_version": 1, "repeat_dictionary": "repeats.json"}, "base.json"
    )
    base = base_settings(base_file)
    assert base["repeat_dictionary"] == str((folder / "repeats.json").resolve())
    (point,) = build_points({"hybrid.seed": [SEED + 1]}, base, "hybrid")
    written = _write(tmp_path, point.config, "overlay.json")
    loaded = load_settings(written)
    assert loaded.hybrid.seed == SEED + 1
    assert loaded.repeat_dictionary == base["repeat_dictionary"]
    assert point.config["schema_version"] == 1


def test_base_settings_defaults_and_invalid_base(tmp_path: Path) -> None:
    assert base_settings(None) == settings_as_dict(DEFAULT_SETTINGS)
    bad = _write(tmp_path, {"schema_version": 1, "run": {"threads": 0}}, "bad.json")
    with pytest.raises(ValueError, match=r"run\.threads"):
        base_settings(bad)


def test_loader_type_errors_become_value_errors(monkeypatch: pytest.MonkeyPatch) -> None:
    import muc_one_span.benchsim.calibration_grid as grid_module

    def broken(path: object) -> None:
        raise TypeError("unexpected keyword")

    monkeypatch.setattr(grid_module, "load_settings", broken)
    with pytest.raises(ValueError, match="unexpected keyword"):
        base_settings(None)
    with pytest.raises(ValueError, match="unexpected keyword"):
        build_points({}, settings_as_dict(DEFAULT_SETTINGS), "hybrid")


def test_a_non_section_key_is_rejected() -> None:
    with pytest.raises(ValueError, match="not a settings section"):
        build_points({"schema_version.x": [1]}, settings_as_dict(DEFAULT_SETTINGS), "hybrid")


def test_stages_and_length_stage_keys_are_all_hybrid_fields() -> None:
    assert STAGES == (DEFAULT_STAGE, LENGTHS_STAGE)
    hybrid_fields = set(settings_as_dict(DEFAULT_SETTINGS)["hybrid"])
    assert LENGTH_STAGE_KEYS
    for key in LENGTH_STAGE_KEYS:
        section, name = key.split(".")
        assert section == "hybrid" and name in hybrid_fields


def _settings_attributes_read(module: ModuleType) -> set[str]:
    """Every ``<name>.<attr>`` read in ``module`` where ``<name>`` is a settings
    parameter (`SETTINGS_PARAM_NAMES`): the `HybridSettings` fields that module's
    source actually accesses, found by walking its AST (not by import/introspection,
    since the fields are read as plain attribute access, not enumerated anywhere)."""
    tree = ast.parse(Path(inspect.getfile(module)).read_text())
    return {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id in SETTINGS_PARAM_NAMES
    }


def test_length_stage_keys_matches_every_setting_spans_lengths_and_smear_read() -> None:
    """`LENGTH_STAGE_KEYS` must track `hybrid.spans`/`hybrid.lengths`/`hybrid.smear`'s
    own settings usage exactly: neither a stale entry (a key that no longer affects the
    length model) nor a missing one (a key that does, silently refused for calibration)
    should be able to drift in without this test failing."""
    read = set().union(
        *(_settings_attributes_read(m) for m in (hybrid_spans, hybrid_lengths, hybrid_smear))
    )
    declared = {key.removeprefix("hybrid.") for key in LENGTH_STAGE_KEYS}
    assert read == declared


def test_check_stage_keys_is_a_noop_for_the_full_stage() -> None:
    check_stage_keys({"hybrid.n_poa": [1]}, DEFAULT_STAGE)  # not a length-model key, but allowed


@pytest.mark.parametrize("key", sorted(LENGTH_STAGE_KEYS)[:3])
def test_check_stage_keys_allows_length_model_keys(key: str) -> None:
    check_stage_keys({key: [1]}, LENGTHS_STAGE)


def test_check_stage_keys_refuses_non_length_keys_for_the_lengths_stage() -> None:
    with pytest.raises(ValueError, match=r"hybrid\.n_poa"):
        check_stage_keys({"hybrid.n_poa": [1], "hybrid.min_span_units": [15]}, LENGTHS_STAGE)
