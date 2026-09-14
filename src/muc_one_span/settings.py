"""Immutable, validated runtime choices shared by library callers and JSON configuration."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any


def _integer(name: str, value: object, minimum: int = 0) -> None:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _number(name: str, value: object, minimum: float = 0, maximum: float | None = None) -> None:
    try:
        finite = isinstance(value, (int, float)) and math.isfinite(value)
    except OverflowError:
        finite = False
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not finite
        or value < minimum
        or (maximum is not None and value > maximum)
    ):
        limit = f" in [{minimum}, {maximum}]" if maximum is not None else f" >= {minimum}"
        raise ValueError(f"{name} must be a finite number{limit}")


def _boolean(name: str, value: object) -> None:
    if type(value) is not bool:
        raise ValueError(f"{name} must be a boolean")


def _string(name: str, value: object, *, optional: bool = False, empty: bool = False) -> None:
    if optional and value is None:
        return
    if not isinstance(value, str) or (not value.strip() and not (empty and value == "")):
        raise ValueError(f"{name} must be a nonempty string" + (" or null" if optional else ""))


def _choice(name: str, value: object, allowed: tuple[str, ...]) -> None:
    if not isinstance(value, str) or value not in allowed:
        raise ValueError(f"{name} must be one of {allowed!r}")


@dataclass(frozen=True)
class RunSettings:
    """Execution defaults; an empty model lets the existing tool choose its default."""

    threads: int = 4
    platform: str = "hifi"
    min_coverage: int = 10
    min_qual: float = 5
    clair3_model: str = ""
    reference: str | None = None
    minimap2_preset: str | None = None
    report: bool = False

    def __post_init__(self) -> None:
        _integer("run.threads", self.threads, 1)
        _integer("run.min_coverage", self.min_coverage, 1)
        _number("run.min_qual", self.min_qual)
        if self.platform not in ("hifi", "ont"):
            raise ValueError("run.platform must be hifi or ont")
        _string("run.clair3_model", self.clair3_model, empty=True)
        _string("run.reference", self.reference, optional=True)
        _string("run.minimap2_preset", self.minimap2_preset, optional=True)
        _boolean("run.report", self.report)


@dataclass(frozen=True)
class AlleleSelectionSettings:
    """Length-distribution separation and refinement heuristics."""

    min_gap: int = 5
    valley_min_points: int = 3
    valley_min_separation: int = 3
    refinement_max_shift: int = 1
    refinement_metric: str = "auto"
    score_margin_hifi: int = 43
    score_margin_ont: int = 42
    min_dominant_reads_hifi: int = 3
    min_dominant_reads_ont: int = 4
    min_dominance_ratio: float = 0.01

    def __post_init__(self) -> None:
        _integer("allele_selection.min_gap", self.min_gap, 1)
        _integer("allele_selection.valley_min_points", self.valley_min_points, 3)
        _integer("allele_selection.valley_min_separation", self.valley_min_separation, 1)
        _integer("allele_selection.refinement_max_shift", self.refinement_max_shift)
        _choice(
            "allele_selection.refinement_metric", self.refinement_metric, ("auto", "as", "indel")
        )
        _integer("allele_selection.score_margin_hifi", self.score_margin_hifi, 0)
        _integer("allele_selection.score_margin_ont", self.score_margin_ont, 0)
        _integer("allele_selection.min_dominant_reads_hifi", self.min_dominant_reads_hifi, 1)
        _integer("allele_selection.min_dominant_reads_ont", self.min_dominant_reads_ont, 1)
        _number("allele_selection.min_dominance_ratio", self.min_dominance_ratio, 0.0, 1.0)


@dataclass(frozen=True)
class ClassificationSettings:
    """Repeat-fitting search limits, segmentation and novelty thresholds."""

    max_indel_probe: int = 30
    max_fit_edit_distance: int = 3
    novel_repeat_edit_distance: int = 2
    minimum_unit_fraction: float = 0.5
    early_stop_edit_distance: int = 1
    strict_segmentation: bool = False

    def __post_init__(self) -> None:
        for name in (
            "max_indel_probe",
            "max_fit_edit_distance",
            "novel_repeat_edit_distance",
            "early_stop_edit_distance",
        ):
            _integer(f"classification.{name}", getattr(self, name))
        _number("classification.minimum_unit_fraction", self.minimum_unit_fraction, 0, 1)
        if self.minimum_unit_fraction == 0:
            raise ValueError("classification.minimum_unit_fraction must be > 0")
        _boolean("classification.strict_segmentation", self.strict_segmentation)


@dataclass(frozen=True)
class ConsensusSettings:
    """Flanking reference extent and exact boundary-anchor search parameters."""

    flank_length: int = 500
    anchor_bases: int = 20
    anchor_tolerance: int = 50
    proximal_flank: bool = True
    haploid_majority: bool = True
    haploid_min_qual: float = 4.0

    def __post_init__(self) -> None:
        _integer("consensus.flank_length", self.flank_length)
        _integer("consensus.anchor_bases", self.anchor_bases, 1)
        _integer("consensus.anchor_tolerance", self.anchor_tolerance)
        _boolean("consensus.proximal_flank", self.proximal_flank)
        _boolean("consensus.haploid_majority", self.haploid_majority)
        _number("consensus.haploid_min_qual", self.haploid_min_qual, 0.0)

    def validate_flanks(self, left: str, right: str, *, flank_length: int | None = None) -> None:
        """Require the effective extent to fit both known dictionary flanks.

        An explicit length overrides this setting; zero permits empty flanks.
        Reject truncation so reference generation and trimming share coordinates.
        """
        extent = self.flank_length if flank_length is None else flank_length
        _integer("consensus.flank_length", extent)
        if extent > len(left) or extent > len(right):
            raise ValueError(
                f"consensus.flank_length={extent} exceeds available dictionary flank "
                f"lengths (left={len(left)}, right={len(right)})"
            )


@dataclass(frozen=True)
class ConfidenceSettings:
    """Quality-weighted confidence and terminal-repeat penalties."""

    qual_low: float = 5
    qual_high: float = 20
    weight_below: float = 0.3
    weight_low: float = 0.5
    weight_high: float = 1
    absent_weight: float = 0.3
    boundary_repeats: int = 3
    boundary_penalty: float = 0.5

    def __post_init__(self) -> None:
        _number("confidence.qual_low", self.qual_low)
        _number("confidence.qual_high", self.qual_high)
        if self.qual_high <= self.qual_low:
            raise ValueError("confidence.qual_high must be greater than qual_low")
        for name in (
            "weight_below",
            "weight_low",
            "weight_high",
            "absent_weight",
            "boundary_penalty",
        ):
            _number(f"confidence.{name}", getattr(self, name), 0, 1)
        _integer("confidence.boundary_repeats", self.boundary_repeats)


@dataclass(frozen=True)
class CallingSettings:
    """VCF sample naming and explicitly experimental read phasing."""

    sample_name: str = "sample"
    read_phase: bool = False
    haploid_majority: bool = True
    haploid_min_qual: float = 4.0

    def __post_init__(self) -> None:
        _string("calling.sample_name", self.sample_name)
        if any(char.isspace() or not char.isprintable() for char in self.sample_name):
            raise ValueError(
                "calling.sample_name must not contain whitespace or control characters"
            )
        _boolean("calling.read_phase", self.read_phase)
        _boolean("calling.haploid_majority", self.haploid_majority)
        _number("calling.haploid_min_qual", self.haploid_min_qual, 0.0)


@dataclass(frozen=True)
class ReadPhasingSettings:
    """Optional WhatsHap overrides; null preserves installed-tool defaults."""

    internal_downsampling: int | None = None
    mapping_quality: int | None = None

    def __post_init__(self) -> None:
        if self.internal_downsampling is not None:
            _integer("read_phasing.internal_downsampling", self.internal_downsampling, 1)
        if self.mapping_quality is not None:
            _integer("read_phasing.mapping_quality", self.mapping_quality)


@dataclass(frozen=True)
class ReferenceLayoutSettings:
    """Ordered fixed repeat IDs, independent of dictionary categories containing alternatives."""

    pre: tuple[str, ...] = ("1", "2", "3", "4", "5")
    after: tuple[str, ...] = ("6", "7", "8", "9")
    min_units: int = 1
    max_units: int = 150

    def __post_init__(self) -> None:
        _integer("reference_layout.min_units", self.min_units, 1)
        _integer("reference_layout.max_units", self.max_units, self.min_units)
        for name, ids in (("pre", self.pre), ("after", self.after)):
            if not isinstance(ids, tuple) or not ids:
                raise ValueError(f"reference_layout.{name} must be a nonempty tuple of IDs")
            for value in ids:
                _string(f"reference_layout.{name} ID", value)
            if len(ids) != len(set(ids)):
                raise ValueError(f"reference_layout.{name} contains duplicate IDs")
        if set(self.pre) & set(self.after):
            raise ValueError("reference_layout pre and after IDs must be disjoint")

    @property
    def fixed_repeat_count(self) -> int:
        """Count only the selected ordered fixed repeats."""
        return len(self.pre) + len(self.after)

    @property
    def left_anchor_id(self) -> str:
        """Return the first repeat next to the left external flank."""
        return self.pre[0]

    @property
    def right_anchor_id(self) -> str:
        """Return the last repeat next to the right external flank."""
        return self.after[-1]

    def validate_repeats(self, repeats: Mapping[str, str]) -> None:
        """Reject a layout whose selected IDs are absent from the loaded dictionary."""
        missing = [name for name in self.pre + self.after if name not in repeats]
        if missing:
            raise ValueError(
                f"reference_layout IDs missing from repeat dictionary: {', '.join(missing)}"
            )


@dataclass(frozen=True)
class RuntimeSettings:
    """Complete schema-one settings; sections remain immutable when passed to workers."""

    schema_version: int = 1
    run: RunSettings = field(default_factory=RunSettings)
    allele_selection: AlleleSelectionSettings = field(default_factory=AlleleSelectionSettings)
    classification: ClassificationSettings = field(default_factory=ClassificationSettings)
    consensus: ConsensusSettings = field(default_factory=ConsensusSettings)
    confidence: ConfidenceSettings = field(default_factory=ConfidenceSettings)
    calling: CallingSettings = field(default_factory=CallingSettings)
    read_phasing: ReadPhasingSettings = field(default_factory=ReadPhasingSettings)
    reference_layout: ReferenceLayoutSettings = field(default_factory=ReferenceLayoutSettings)
    repeat_dictionary: str | None = None

    def __post_init__(self) -> None:
        if type(self.schema_version) is not int or self.schema_version != 1:
            raise ValueError("schema_version must be integer 1")
        for name, expected in _SECTIONS.items():
            if not isinstance(getattr(self, name), expected):
                raise ValueError(f"{name} must be {expected.__name__}")
        _string("repeat_dictionary", self.repeat_dictionary, optional=True)


_SECTIONS = {
    "run": RunSettings,
    "allele_selection": AlleleSelectionSettings,
    "classification": ClassificationSettings,
    "consensus": ConsensusSettings,
    "confidence": ConfidenceSettings,
    "calling": CallingSettings,
    "read_phasing": ReadPhasingSettings,
    "reference_layout": ReferenceLayoutSettings,
}
DEFAULT_SETTINGS = RuntimeSettings()
DEFAULT_LAYOUT = DEFAULT_SETTINGS.reference_layout


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate configuration key: {key}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise ValueError(f"Nonfinite JSON number: {value}")


def _resolve_path(value: str | None, folder: Path) -> str | None:
    if value is None or value == "":
        return value
    path = Path(value)
    return str((folder / path).resolve())


def load_settings(path: Path | None) -> RuntimeSettings:
    """Load strict schema-one JSON, resolving resource paths relative to its location.

    Missing sections/fields use central defaults. A supplied file must declare
    schema_version=1; duplicate/unknown keys and invalid values raise ValueError.
    Filesystem read failures propagate. No resources or external tools are loaded.
    """
    if path is None:
        return DEFAULT_SETTINGS
    data = json.loads(
        path.read_text(), object_pairs_hook=_unique_object, parse_constant=_reject_constant
    )
    if not isinstance(data, dict):
        raise ValueError("Runtime configuration must be a JSON object")
    unknown = data.keys() - {f.name for f in fields(RuntimeSettings)}
    if unknown:
        raise ValueError(f"Unknown configuration fields: {', '.join(sorted(unknown))}")
    if "schema_version" not in data:
        raise ValueError("Runtime configuration requires schema_version=1")
    for name, constructor in _SECTIONS.items():
        if name not in data:
            continue
        values = data[name]
        if not isinstance(values, dict):
            raise ValueError(f"{name} must be a JSON object")
        unknown = values.keys() - {f.name for f in fields(constructor)}
        if unknown:
            raise ValueError(f"Unknown {name} fields: {', '.join(sorted(unknown))}")
        if name == "reference_layout":
            for key in ("pre", "after"):
                if key not in values:
                    continue
                value = values[key]
                if not isinstance(value, list):
                    raise ValueError(f"reference_layout.{key} must be a JSON array")
                values[key] = tuple(value)
        data[name] = constructor(**values)
    settings = RuntimeSettings(**data)
    # Validate raw string types before converting paths, so coercion never hides errors.
    run = asdict(settings.run)
    run["reference"] = _resolve_path(settings.run.reference, path.resolve().parent)
    run["clair3_model"] = _resolve_path(settings.run.clair3_model, path.resolve().parent)
    return RuntimeSettings(
        **{
            **{f.name: getattr(settings, f.name) for f in fields(RuntimeSettings)},
            "run": RunSettings(**run),
            "repeat_dictionary": _resolve_path(settings.repeat_dictionary, path.resolve().parent),
        }
    )


def settings_as_dict(settings: RuntimeSettings) -> dict[str, Any]:
    """Return a JSON-compatible effective configuration suitable for provenance/roundtrip."""
    result = asdict(settings)
    result["reference_layout"]["pre"] = list(settings.reference_layout.pre)
    result["reference_layout"]["after"] = list(settings.reference_layout.after)
    return result
