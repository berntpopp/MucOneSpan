"""Immutable, validated runtime choices shared by library callers and JSON configuration."""

from __future__ import annotations

import json
import logging
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

from muc_one_span.settings_hybrid import SMEAR_CORRECTIONS, HybridSettings
from muc_one_span.settings_validation import (
    _boolean,
    _choice,
    _integer,
    _number,
    _open_unit_interval,
    _string,
)

# Re-exported for backward compatibility: these validation helpers and HybridSettings
# used to live in this module; existing `from muc_one_span.settings import ...` call
# sites (in src and tests) keep working unchanged. See settings_validation.py and
# settings_hybrid.py for their definitions.
__all__ = [
    "SMEAR_CORRECTIONS",
    "HybridSettings",
    "_boolean",
    "_choice",
    "_integer",
    "_number",
    "_open_unit_interval",
    "_string",
]

logger = logging.getLogger(__name__)


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
    mapping_timeout: float = 3600.0
    report: bool = False
    report_igv: str = "off"
    engine: str = "ladder"
    assay: str = "amplicon"

    def __post_init__(self) -> None:
        _integer("run.threads", self.threads, 1)
        _integer("run.min_coverage", self.min_coverage, 1)
        _number("run.min_qual", self.min_qual)
        _number("run.mapping_timeout", self.mapping_timeout)
        if self.mapping_timeout == 0:
            raise ValueError("run.mapping_timeout must be > 0")
        if self.platform not in ("hifi", "ont"):
            raise ValueError("run.platform must be hifi or ont")
        _string("run.clair3_model", self.clair3_model, empty=True)
        _string("run.reference", self.reference, optional=True)
        _string("run.minimap2_preset", self.minimap2_preset, optional=True)
        _boolean("run.report", self.report)
        if self.report_igv not in ("off", "embedded", "sidecar"):
            raise ValueError("run.report_igv must be off, embedded, or sidecar")
        _choice("run.engine", self.engine, ("ladder", "hybrid"))
        _choice("run.assay", self.assay, ("amplicon", "genomic"))


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
    secondary_mode_min_fraction: float = 0.2
    min_allele_primary_records: int = 30
    # Ladder heuristics formerly hardcoded (#74); defaults reproduce v0.16.0.
    refinement_min_supported_records: int = 3
    refinement_supported_fraction: float = 0.25
    refinement_min_shift_ont: int = 2
    valley_min_canonical_repeats: int = 10
    minority_min_alignment_records: int = 3
    dominance_close_candidate_repeats: int = 6
    dominance_zero_primary_extra_reads: int = 2
    read_length_split_min_reads: int = 5
    read_length_split_min_fraction: float = 0.15
    read_length_split_bin_bp: int = 5
    read_length_split_min_delta_bp: int = 45
    read_length_split_max_delta_bp: int = 320
    read_length_split_unit_tolerance_bp: int = 15
    read_length_split_offset_bp: int = 30

    def __post_init__(self) -> None:
        self._validate_ladder_heuristics()
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
        _number(
            "allele_selection.secondary_mode_min_fraction",
            self.secondary_mode_min_fraction,
            0.0,
            1.0,
        )
        if self.secondary_mode_min_fraction == 0:
            raise ValueError("allele_selection.secondary_mode_min_fraction must be > 0")
        _integer("allele_selection.min_allele_primary_records", self.min_allele_primary_records, 1)

    def _validate_ladder_heuristics(self) -> None:
        """Integers >= their minimum, fractions in (0, 1], and an ordered delta window."""
        for name, minimum in _LADDER_INTEGER_MINIMUMS:
            _integer(f"allele_selection.{name}", getattr(self, name), minimum)
        for name in ("refinement_supported_fraction", "read_length_split_min_fraction"):
            _number(f"allele_selection.{name}", getattr(self, name), 0.0, 1.0)
            if getattr(self, name) == 0:
                raise ValueError(f"allele_selection.{name} must be > 0")
        if self.read_length_split_max_delta_bp <= self.read_length_split_min_delta_bp:
            raise ValueError(
                "allele_selection.read_length_split_max_delta_bp must be greater than "
                "read_length_split_min_delta_bp"
            )


_LADDER_INTEGER_MINIMUMS = (
    ("refinement_min_supported_records", 1),
    ("refinement_min_shift_ont", 0),
    ("valley_min_canonical_repeats", 1),
    ("minority_min_alignment_records", 1),
    ("dominance_close_candidate_repeats", 1),
    ("dominance_zero_primary_extra_reads", 1),
    ("read_length_split_min_reads", 1),
    ("read_length_split_bin_bp", 1),
    ("read_length_split_min_delta_bp", 1),
    ("read_length_split_max_delta_bp", 1),
    ("read_length_split_unit_tolerance_bp", 0),
    ("read_length_split_offset_bp", 0),
)


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
    """Flanking reference extent and exact boundary-anchor search parameters.

    ``haploid_majority`` and ``haploid_min_qual`` are deprecated no-ops kept for
    configuration compatibility; ``calling.haploid_*`` controls haploid calling.
    """

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
    haploid_min_qual: float | None = 4.0
    haploid_alt_fraction: float = 0.5
    haploid_ref_fraction: float = 0.2
    stage_discordance_min_af: float = 0.5
    stage_discordance_min_depth: int = 10

    def __post_init__(self) -> None:
        _string("calling.sample_name", self.sample_name)
        if any(char.isspace() or not char.isprintable() for char in self.sample_name):
            raise ValueError(
                "calling.sample_name must not contain whitespace or control characters"
            )
        _boolean("calling.read_phase", self.read_phase)
        _boolean("calling.haploid_majority", self.haploid_majority)
        if self.haploid_min_qual is not None:
            _number("calling.haploid_min_qual", self.haploid_min_qual, 0.0)
        _number("calling.haploid_alt_fraction", self.haploid_alt_fraction, 0.0, 1.0)
        _number("calling.haploid_ref_fraction", self.haploid_ref_fraction, 0.0, 1.0)
        if self.haploid_ref_fraction >= self.haploid_alt_fraction:
            raise ValueError(
                "calling.haploid_ref_fraction must be below calling.haploid_alt_fraction"
            )
        _number("calling.stage_discordance_min_af", self.stage_discordance_min_af, 0.0, 1.0)
        if self.stage_discordance_min_af == 0:
            raise ValueError("calling.stage_discordance_min_af must be > 0")
        _integer("calling.stage_discordance_min_depth", self.stage_discordance_min_depth, 1)


@dataclass(frozen=True)
class ReadPhasingSettings:
    """Optional WhatsHap overrides; null preserves installed-tool defaults."""

    internal_downsampling: int | None = None
    mapping_quality: int | None = None
    # Haplotag split: minimum reads per haplotype (formerly ``min_dp``, fixed at 5).
    min_haplotype_reads: int = 5

    def __post_init__(self) -> None:
        _integer("read_phasing.min_haplotype_reads", self.min_haplotype_reads, 1)
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
class ClinicalDecisionSettings:
    """Thresholds for the report's final clinical decision banner (``report.py``).

    Defaults reproduce v0.16.0 exactly: a summed ambiguous-base count above
    ``max_ambiguous_bases`` and, for legacy summaries without per-allele depth,
    a total read count below ``legacy_min_total_reads`` each block NEGATIVE.
    """

    max_ambiguous_bases: int = 10
    legacy_min_total_reads: int = 30

    def __post_init__(self) -> None:
        _integer("clinical_decision.max_ambiguous_bases", self.max_ambiguous_bases, 0)
        _integer("clinical_decision.legacy_min_total_reads", self.legacy_min_total_reads, 1)


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
    hybrid: HybridSettings = field(default_factory=HybridSettings)
    clinical_decision: ClinicalDecisionSettings = field(default_factory=ClinicalDecisionSettings)
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
    "hybrid": HybridSettings,
    "clinical_decision": ClinicalDecisionSettings,
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


def _coerce_reference_layout_lists(values: dict[str, Any]) -> None:
    """Turn JSON arrays into the hashable tuples ``ReferenceLayoutSettings`` requires."""
    for key in ("pre", "after"):
        if key not in values:
            continue
        value = values[key]
        if not isinstance(value, list):
            raise ValueError(f"reference_layout.{key} must be a JSON array")
        values[key] = tuple(value)


def build_settings_section(
    name: str,
    constructor: type,
    values: object,
    *,
    preprocess: Callable[[dict[str, Any]], None] | None = None,
) -> Any:
    """Build one settings section from a JSON object: shared by a configuration file's
    sections (``load_settings``) and a recorded ``configuration.settings`` section
    (``decision_settings.resolve_decision_settings``).

    A field missing from *values* keeps the dataclass default. An unrecognised
    field always raises ``ValueError`` (never a bare ``TypeError`` from unpacking
    a non-mapping, and never silently ignored or defaulted).
    """
    if not isinstance(values, dict):
        raise ValueError(f"{name} must be a JSON object")
    unknown = values.keys() - {f.name for f in fields(constructor)}
    if unknown:
        raise ValueError(f"Unknown {name} fields: {', '.join(sorted(unknown))}")
    if preprocess is not None:
        preprocess(values)
    return constructor(**values)


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
        preprocess = _coerce_reference_layout_lists if name == "reference_layout" else None
        data[name] = build_settings_section(name, constructor, data[name], preprocess=preprocess)
    settings = RuntimeSettings(**data)
    consensus, legacy = settings.consensus, ConsensusSettings()
    if (consensus.haploid_majority, consensus.haploid_min_qual) != (
        legacy.haploid_majority,
        legacy.haploid_min_qual,
    ):
        # Logged rather than a DeprecationWarning, which default filters hide from CLI users.
        logger.warning(
            "Deprecated: consensus.haploid_majority and consensus.haploid_min_qual have no "
            "effect; use calling.haploid_majority and calling.haploid_min_qual."
        )
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
