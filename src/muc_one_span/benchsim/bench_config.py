"""Validated MucSim-Bench settings: every tunable number of the benchmark in one place.

The dataclass defaults are the documented defaults (see ``docs/benchmark.md``,
"Configuration"). `load_bench_config` overlays a strict JSON file on them:
unknown sections or fields, duplicate keys and invalid values raise
``ValueError``. The effective settings are hashed (`BenchConfig.sha256`) and
recorded with generated cases and reports.

Domain constants are not settings: the repeat-unit length and the conserved
repeat IDs come from the bundled repeat dictionary, and the conserved head and
tail positions from `muc_one_span.settings.DEFAULT_LAYOUT` (see `design`).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

from muc_one_span.settings import DEFAULT_SETTINGS

from .bench_checks import ERROR_LEVEL_NAMES, PCR_LEVEL_NAMES
from .bench_checks import check_int as _int
from .bench_checks import check_names as _names
from .bench_checks import check_num as _num
from .bench_sets import SetsConfig

__all__ = [
    "ERROR_LEVEL_NAMES",
    "PCR_LEVEL_NAMES",
    "TARGET_BASIS_NAMES",
    "TARGET_COMPARATOR_NAMES",
    "TARGET_METRIC_NAMES",
    "SetsConfig",
    "Target",
    "TargetsConfig",
]

SCHEMA_VERSION = 1
PERCENT = 100  # probability -> percent (unit conversion)
# Names with semantics implemented elsewhere (`structures`, `generate`, `design`).
COMPOSITION_NAMES = ("markov", "real_derived", "rare_units")
DELTA_CLASS_NAMES = ("0_identical", "0_different", "1", "2", "3-5", "6-20", ">20")
EQUAL_LENGTH_CLASSES = ("0_identical", "0_different")  # length difference must be 0
# Clinical decision states (`muc_one_span.report.compute_clinical_decision`, plus the
# evaluator's NO_CALL for unreadable or failed runs).
DECISION_NAMES = ("PATHOGENIC", "INCONCLUSIVE", "NO_PATHOGENIC_VARIANT_DETECTED", "NO_CALL")
# Design factors the reason atlas can stratify by (`report.normalize_rows` row keys).
ATLAS_STRATUM_NAMES = (
    "profile",
    "event",
    "event_position",
    "delta_class",
    "depth",
    "composition",
    "pcr",
    "smear",
    "chimera",
    "error",
)
# Depth compared with the atlas depth gate: the design target, or the lowest
# realized spanning depth over the case's alleles (``case.json`` realized_depth).
DEPTH_BASIS_NAMES = ("design", "realized_min_allele")
# Owner-approved absolute targets (task 12e, `targets` section): metrics a bench set's
# `targets.by_set` entry may name (`benchsim.targets`), comparators ("ge" >=, "le" <=,
# candidate value vs. threshold) and how a target is judged (`benchsim.targets.evaluate_targets`).
TARGET_METRIC_NAMES = ("pathogenic_rate", "inconclusive_rate", "false_positive_rate")
TARGET_COMPARATOR_NAMES = ("ge", "le")
TARGET_BASIS_NAMES = ("point", "ci_bound")
# Sections that shape generated cases (designs, amounts, read profiles, structures).
# With the case's own set levels (`BenchConfig.generation_sha256`) their hash decides
# whether `generate` may reuse a case; report, realism, run and atlas settings, set
# names, descriptions and other sets' or profiles' levels do not.
GENERATION_SECTIONS = ("design", "amount", "profiles", "structures")
# Run-time knobs inside generation sections, excluded from the generation hash.
RUNTIME_ONLY = {"profiles": ("simulator_threads",)}
_WEIGHT_SUM_TOL = 1e-9  # float round-off allowed when composition weights sum to 1


@dataclass(frozen=True)
class DesignConfig:
    """Split sizes and the biological design factors (spec section 5).

    Technical factor levels (depth, PCR, error, smear, chimera) are per benchmark
    set and profile (`bench_sets.SetsConfig`).
    """

    split_sizes: dict[str, int] = field(
        default_factory=lambda: {"dev": 300, "val": 300, "test": 800, "stress": 100}
    )
    normal_fraction: float = 0.35
    length_min: int = 20
    length_max: int = 130
    delta_ranges: dict[str, tuple[int, int]] = field(
        default_factory=lambda: {
            "0_identical": (0, 0),
            "0_different": (0, 0),
            "1": (1, 1),
            "2": (2, 2),
            "3-5": (3, 5),
            "6-20": (6, 20),
            ">20": (21, 90),
        }
    )
    compositions: dict[str, float] = field(
        default_factory=lambda: {"markov": 0.75, "real_derived": 0.20, "rare_units": 0.05}
    )
    # "first10"/"last10": the event lies in this leading/trailing fraction of the allele.
    position_fraction: float = 0.1

    def __post_init__(self) -> None:
        if not self.split_sizes:
            raise ValueError("design.split_sizes must name at least one split")
        for split, n in self.split_sizes.items():
            _int(f"design.split_sizes.{split}", n, 1)
        _num("design.normal_fraction", self.normal_fraction, 0, 1)
        _int("design.length_min", self.length_min, 1)
        _int("design.length_max", self.length_max, self.length_min)
        self._check_delta_ranges()
        self._check_compositions()
        _num("design.position_fraction", self.position_fraction, 0, 0.5, open_lo=True)

    def _check_delta_ranges(self) -> None:
        if not self.delta_ranges:
            raise ValueError("design.delta_ranges must name at least one class")
        span = self.length_max - self.length_min
        for name, bounds in self.delta_ranges.items():
            key = f"design.delta_ranges.{name}"
            if name not in DELTA_CLASS_NAMES:
                raise ValueError(f"{key}: unknown class (known: {', '.join(DELTA_CLASS_NAMES)})")
            if not isinstance(bounds, tuple) or len(bounds) != 2:
                raise ValueError(f"{key} must be a [min, max] pair")
            lo, hi = bounds
            _int(key, lo, 0)
            _int(key, hi, lo)
            if hi > span:
                raise ValueError(f"{key} exceeds the length range")
            if name in EQUAL_LENGTH_CLASSES and hi != 0:
                raise ValueError(f"{key} must be [0, 0] (equal allele lengths)")

    def _check_compositions(self) -> None:
        if not self.compositions:
            raise ValueError("design.compositions must name at least one composition")
        for name, weight in self.compositions.items():
            key = f"design.compositions.{name}"
            if name not in COMPOSITION_NAMES:
                raise ValueError(
                    f"{key}: unknown composition (known: {', '.join(COMPOSITION_NAMES)})"
                )
            _num(key, weight, 0, 1)
        weights = self.compositions.values()
        if not any(w > 0 for w in weights):
            raise ValueError("design.compositions needs at least one positive weight")
        if abs(sum(weights) - 1) > _WEIGHT_SUM_TOL:
            raise ValueError("design.compositions weights must sum to 1")


@dataclass(frozen=True)
class AmountConfig:
    """Template and read amounts that reach the design depth (`depth`, `generate`)."""

    # PCR share of the minor allele: 1 / (1 + exp(slope * length difference in units)).
    pcr_slope_per_unit: dict[str, float] = field(
        default_factory=lambda: {"calibrated": 0.056, "strong": 0.112, "none": 0.0}
    )
    # Floor on the minor-allele share used to size amplicon templates; below it the
    # template count explodes and the minor allele gets below-target depth instead.
    min_minor_share: float = 0.05
    genomic_mc_draws: int = 20000
    # MucOneUp FragmentModel defaults, used when a base profile has no fragments block.
    fragment_length_median: float = 5000.0
    fragment_length_sigma: float = 0.5

    def __post_init__(self) -> None:
        if set(self.pcr_slope_per_unit) != set(PCR_LEVEL_NAMES):
            raise ValueError(f"amount.pcr_slope_per_unit must define {PCR_LEVEL_NAMES!r}")
        for name, slope in self.pcr_slope_per_unit.items():
            _num(f"amount.pcr_slope_per_unit.{name}", slope, 0)
        _num("amount.min_minor_share", self.min_minor_share, 0, 0.5, open_lo=True)
        _int("amount.genomic_mc_draws", self.genomic_mc_draws, 1)
        _num("amount.fragment_length_median", self.fragment_length_median, 1)
        _num("amount.fragment_length_sigma", self.fragment_length_sigma, 0, open_lo=True)


@dataclass(frozen=True)
class ProfileConfig:
    """Error and PCR levels layered on the MucOneUp base profiles (`profiles`)."""

    r10_pcr_alpha: float = 9.27e-5  # madritsch2025_r10 preset alpha
    strong_pcr_alpha_factor: float = 2.0
    poor_error_scale: float = 1.5
    hifi_poor_accuracy_mean: float = 0.95
    # Threads per case for MucOneUp's simulator tools (pbsim3, ccs, samtools,
    # minimap2), written into every profile variant; generate --jobs times this
    # is the approximate core use of a generate run.
    simulator_threads: int = 2

    def __post_init__(self) -> None:
        _num("profiles.r10_pcr_alpha", self.r10_pcr_alpha, 0, open_lo=True)
        _num("profiles.strong_pcr_alpha_factor", self.strong_pcr_alpha_factor, 1)
        _num("profiles.poor_error_scale", self.poor_error_scale, 1)
        _num("profiles.hifi_poor_accuracy_mean", self.hifi_poor_accuracy_mean, 0, 1, open_lo=True)
        _int("profiles.simulator_threads", self.simulator_threads, 1)


@dataclass(frozen=True)
class StructureConfig:
    """Rare-unit structures (`structures`)."""

    rare_fraction: float = 0.10  # share of interior units replaced
    rare_usage_max: float = 0.01  # a unit is rare below this Markov visit frequency
    stationary_steps: int = 2000  # power iterations for the visit frequency

    def __post_init__(self) -> None:
        _num("structures.rare_fraction", self.rare_fraction, 0, 1, open_lo=True)
        _num("structures.rare_usage_max", self.rare_usage_max, 0, 1, open_lo=True)
        _int("structures.stationary_steps", self.stationary_steps, 1)


@dataclass(frozen=True)
class RealismConfig:
    """Realism metric definitions and pass tolerances (spec section 4)."""

    error_rel_tol: float = 0.20
    c7_abs_tol: float = 0.03
    range_median_rel_tol: float = 0.25
    jsd_max: float = 0.1
    slope_abs_tol: float = 0.01
    # Span-offset histogram; must match the target file's bin_lo_bp (checked in compare).
    offset_bin_bp: int = 15
    offset_bin_min: int = -12
    offset_bin_max: int = 4
    size_split_units: int = 55  # histogram keys lt<N>u / ge<N>u
    c7_run_length: int = 7
    c7_extend_max: int = 3  # adjacent read C counted beyond the aligned run
    on_peak_min_bp: float = 30.0
    on_peak_rel: float = 0.012
    off_peak_units: float = 1.5

    def __post_init__(self) -> None:
        for name in ("error_rel_tol", "c7_abs_tol", "range_median_rel_tol", "slope_abs_tol"):
            _num(f"realism.{name}", getattr(self, name), 0)
        _num("realism.jsd_max", self.jsd_max, 0, 1)
        _int("realism.offset_bin_bp", self.offset_bin_bp, 1)
        if type(self.offset_bin_min) is not int or type(self.offset_bin_max) is not int:
            raise ValueError("realism.offset_bin_min/max must be integers")
        if self.offset_bin_min >= self.offset_bin_max:
            raise ValueError("realism.offset_bin_min must be below offset_bin_max")
        _int("realism.size_split_units", self.size_split_units, 1)
        _int("realism.c7_run_length", self.c7_run_length, 1)
        _int("realism.c7_extend_max", self.c7_extend_max, 0)
        _num("realism.on_peak_min_bp", self.on_peak_min_bp, 0)
        _num("realism.on_peak_rel", self.on_peak_rel, 0)
        _num("realism.off_peak_units", self.off_peak_units, 0, open_lo=True)

    @property
    def n_bins(self) -> int:
        """Number of span-offset histogram bins."""
        return self.offset_bin_max - self.offset_bin_min + 1

    @property
    def size_keys(self) -> tuple[str, str]:
        """Histogram keys below / at-or-above `size_split_units`."""
        return f"lt{self.size_split_units}u", f"ge{self.size_split_units}u"

    def bin_lo_bp(self) -> list[int]:
        """Lower edge (bp) of each bin; bin ``b`` holds ``(d + bp // 2) // bp == b``."""
        half = self.offset_bin_bp // 2
        return [
            b * self.offset_bin_bp - half
            for b in range(self.offset_bin_min, self.offset_bin_max + 1)
        ]


@dataclass(frozen=True)
class ReportConfig:
    """Decision-rule and interval parameters (spec section 6, `report`)."""

    alpha: float = 0.05
    ni_margin: float = 0.005  # non-inferiority margin on the FP rate difference
    bootstrap_replicates: int = 2000
    bootstrap_seed: int = 0

    def __post_init__(self) -> None:
        _num("report.alpha", self.alpha, 0, 0.5, open_lo=True)
        _num("report.ni_margin", self.ni_margin, 0, 1, open_lo=True)
        _int("report.bootstrap_replicates", self.bootstrap_replicates, 1)
        _int("report.bootstrap_seed", self.bootstrap_seed, 0)


@dataclass(frozen=True)
class Target:
    """One absolute pass/fail check (task 12e): the judged value `comparator` `threshold`.

    ``comparator`` is ``"ge"`` (>=, a rate floor such as PATHOGENIC-on-pathogenic) or
    ``"le"`` (<=, a rate ceiling such as INCONCLUSIVE or false-positive). Validated by
    `_check_targets`, which has the enclosing bench set and metric name for its errors.
    """

    comparator: str
    threshold: float


def _check_targets(set_name: str, metrics: Any) -> None:
    key = f"targets.by_set.{set_name}"
    if not isinstance(metrics, dict):
        raise ValueError(f"{key} must be a JSON object")
    for metric, target in metrics.items():
        where = f"{key}.{metric}"
        if metric not in TARGET_METRIC_NAMES:
            raise ValueError(f"{where}: unknown metric (known: {', '.join(TARGET_METRIC_NAMES)})")
        if not isinstance(target, Target):
            raise ValueError(f"{where} must be a target object")
        if target.comparator not in TARGET_COMPARATOR_NAMES:
            raise ValueError(f"{where}.comparator must be one of {TARGET_COMPARATOR_NAMES!r}")
        _num(f"{where}.threshold", target.threshold, 0, 1)


@dataclass(frozen=True)
class TargetsConfig:
    """Owner-approved absolute targets per bench set (task 12e, spec section 6, `targets`).

    ``basis`` decides how a target is judged (`benchsim.targets.evaluate_targets`): the
    pooled/per-profile point estimate (``"point"``), or the Clopper-Pearson CI bound on
    the side `comparator` cares about (``"ci_bound"``: the lower bound for ``"ge"``, the
    upper bound for ``"le"``). ``by_set`` names, for each bench set, the metrics that gate
    `report.decide`'s adoption verdict; a set absent from `by_set` (or mapped to an empty
    object) is reported without a target, for example `stress`. Like `sets.definitions`,
    overriding `by_set` in a bench-config file replaces the whole map.
    """

    basis: str = "point"
    by_set: dict[str, dict[str, Target]] = field(
        default_factory=lambda: {
            "standard": {
                "pathogenic_rate": Target("ge", 0.80),
                "inconclusive_rate": Target("le", 0.20),
                "false_positive_rate": Target("le", 0.0),
            },
            "clean": {
                "pathogenic_rate": Target("ge", 0.90),
                "inconclusive_rate": Target("le", 0.10),
                "false_positive_rate": Target("le", 0.0),
            },
        }
    )

    def __post_init__(self) -> None:
        if self.basis not in TARGET_BASIS_NAMES:
            raise ValueError(f"targets.basis must be one of {TARGET_BASIS_NAMES!r}")
        if not isinstance(self.by_set, dict):
            raise ValueError("targets.by_set must be a JSON object")
        for set_name, metrics in self.by_set.items():
            _check_targets(set_name, metrics)


@dataclass(frozen=True)
class AtlasConfig:
    """Reason atlas of non-definitive decisions (`atlas`, docs/benchmark.md).

    A case is *expected* non-definitive when its split is listed in
    ``expected_inconclusive_splits`` or its depth (``depth_basis``) is below
    ``min_resolvable_depth``; every other atlas case counts as resolvable.
    """

    decisions: tuple[str, ...] = ("INCONCLUSIVE",)
    strata: tuple[str, ...] = ("depth", "smear", "chimera", "delta_class", "event_position")
    expected_inconclusive_splits: tuple[str, ...] = ("stress",)
    expected_inconclusive_sets: tuple[str, ...] = ("stress",)
    # Default: the caller's per-allele primary-alignment gate for a negative call.
    min_resolvable_depth: int = DEFAULT_SETTINGS.allele_selection.min_allele_primary_records
    depth_basis: str = "realized_min_allele"
    top_reasons: int = 10  # reasons shown per stratum table in report.md (JSON keeps all)

    def __post_init__(self) -> None:
        _names("atlas.decisions", self.decisions, DECISION_NAMES)
        _names("atlas.strata", self.strata, ATLAS_STRATUM_NAMES)
        for name in ("expected_inconclusive_splits", "expected_inconclusive_sets"):
            values = getattr(self, name)
            if not isinstance(values, tuple) or not all(isinstance(v, str) for v in values):
                raise ValueError(f"atlas.{name} must be a list of names")
            if len(set(values)) != len(values):
                raise ValueError(f"atlas.{name} must be distinct")
        _int("atlas.min_resolvable_depth", self.min_resolvable_depth, 1)
        if self.depth_basis not in DEPTH_BASIS_NAMES:
            raise ValueError(f"atlas.depth_basis must be one of {DEPTH_BASIS_NAMES!r}")
        _int("atlas.top_reasons", self.top_reasons, 1)


@dataclass(frozen=True)
class RunConfig:
    """Engine run defaults (`run`)."""

    threads: int = 4

    def __post_init__(self) -> None:
        _int("run.threads", self.threads, 1)


@dataclass(frozen=True)
class BenchConfig:
    """All MucSim-Bench settings; immutable and picklable for worker processes."""

    schema_version: int = SCHEMA_VERSION
    design: DesignConfig = field(default_factory=DesignConfig)
    amount: AmountConfig = field(default_factory=AmountConfig)
    profiles: ProfileConfig = field(default_factory=ProfileConfig)
    structures: StructureConfig = field(default_factory=StructureConfig)
    realism: RealismConfig = field(default_factory=RealismConfig)
    report: ReportConfig = field(default_factory=ReportConfig)
    run: RunConfig = field(default_factory=RunConfig)
    atlas: AtlasConfig = field(default_factory=AtlasConfig)
    sets: SetsConfig = field(default_factory=SetsConfig)
    targets: TargetsConfig = field(default_factory=TargetsConfig)

    def __post_init__(self) -> None:
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"bench config schema_version must be {SCHEMA_VERSION}")
        unknown = set(self.atlas.expected_inconclusive_splits) - set(self.design.split_sizes)
        if unknown:
            raise ValueError(
                "atlas.expected_inconclusive_splits must name splits in design.split_sizes "
                f"(unknown: {', '.join(sorted(unknown))})"
            )
        unknown = set(self.atlas.expected_inconclusive_sets) - set(self.sets.definitions)
        if unknown:
            raise ValueError(
                "atlas.expected_inconclusive_sets must name sets in sets.definitions "
                f"(unknown: {', '.join(sorted(unknown))})"
            )
        unknown = set(self.targets.by_set) - set(self.sets.definitions)
        if unknown:
            raise ValueError(
                "targets.by_set must name sets in sets.definitions "
                f"(unknown: {', '.join(sorted(unknown))})"
            )

    def to_dict(self) -> dict[str, Any]:
        """JSON-ready settings (tuples become lists)."""
        data: dict[str, Any] = json.loads(json.dumps(asdict(self)))
        return data

    def sha256(self) -> str:
        """SHA-256 of the canonical JSON of the effective settings (provenance)."""
        return _digest(self.to_dict())

    def generation_sha256(self, bench_set: str | None = None, profile: str | None = None) -> str:
        """SHA-256 of what shapes one case (the `generate` reuse check).

        Covers the `GENERATION_SECTIONS` without their `RUNTIME_ONLY` knobs and,
        with ``profile``, the levels of ``bench_set`` for that profile (``None``
        when the set or profile is undefined, as for pre-set designs).
        """
        data = self.to_dict()
        payload = {k: data[k] for k in ("schema_version", *GENERATION_SECTIONS)}
        for section, keys in RUNTIME_ONLY.items():
            payload[section] = {k: v for k, v in payload[section].items() if k not in keys}
        if profile is not None:
            definition = data["sets"]["definitions"].get(bench_set) if bench_set else None
            levels = (definition or {}).get("profiles", {}).get(profile)
            payload["set_levels"] = {"set": bench_set, "profile": profile, "levels": levels}
        return _digest(payload)


def _digest(data: dict[str, Any]) -> str:
    text = json.dumps(data, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()


DEFAULT_BENCH_CONFIG = BenchConfig()
_SECTIONS = {f.name: f for f in fields(BenchConfig) if f.name != "schema_version"}


def _unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in pairs:
        if key in out:
            raise ValueError(f"duplicate bench config key: {key}")
        out[key] = value
    return out


def _coerce(default: Any, value: Any, name: str) -> Any:
    """JSON value shaped like the default: lists -> tuples, dict values likewise.

    A dataclass default (a benchmark set, its per-profile levels) needs a JSON
    object with every field and no others.
    """
    if is_dataclass(default) and not isinstance(default, type):
        if not isinstance(value, dict):
            raise ValueError(f"{name} must be a JSON object")
        known = {f.name for f in fields(default)}
        if value.keys() - known:
            raise ValueError(f"unknown {name} fields: {', '.join(sorted(value.keys() - known))}")
        if known - value.keys():
            raise ValueError(f"{name}: missing fields: {', '.join(sorted(known - value.keys()))}")
        coerced = {k: _coerce(getattr(default, k), value[k], f"{name}.{k}") for k in known}
        return type(default)(**coerced)
    if isinstance(default, tuple):
        if not isinstance(value, list):
            raise ValueError(f"{name} must be a JSON array")
        return tuple(value)
    if isinstance(default, dict):
        if not isinstance(value, dict):
            raise ValueError(f"{name} must be a JSON object")
        sample = next(iter(default.values()), None)
        return {k: _coerce(sample, v, f"{name}.{k}") for k, v in value.items()}
    return value


def load_bench_config(path: Path | None = None) -> BenchConfig:
    """Defaults, or defaults overlaid with a strict JSON file (``schema_version`` 1)."""
    if path is None:
        return DEFAULT_BENCH_CONFIG
    data = json.loads(path.read_text(), object_pairs_hook=_unique)
    if not isinstance(data, dict):
        raise ValueError("bench config must be a JSON object")
    unknown = data.keys() - {f.name for f in fields(BenchConfig)}
    if unknown:
        raise ValueError(f"unknown bench config sections: {', '.join(sorted(unknown))}")
    if data.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"bench config requires schema_version={SCHEMA_VERSION}")
    sections: dict[str, Any] = {}
    for name in _SECTIONS:
        base = getattr(DEFAULT_BENCH_CONFIG, name)
        values = data.get(name, {})
        if not isinstance(values, dict):
            raise ValueError(f"{name} must be a JSON object")
        known = {f.name for f in fields(base)}
        if values.keys() - known:
            raise ValueError(f"unknown {name} fields: {', '.join(sorted(values.keys() - known))}")
        merged = {
            key: _coerce(getattr(base, key), values[key], f"{name}.{key}")
            if key in values
            else getattr(base, key)
            for key in known
        }
        try:
            sections[name] = type(base)(**merged)
        except TypeError as exc:
            raise ValueError(f"{name}: invalid value type ({exc})") from exc
    return BenchConfig(schema_version=SCHEMA_VERSION, **sections)
