"""Immutable, validated runtime choices shared by library callers and JSON configuration."""

from __future__ import annotations

import json
import logging
import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


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


def _open_unit_interval(name: str, value: object) -> None:
    _number(name, value, 0, 1)
    if value in (0, 1):
        raise ValueError(f"{name} must be strictly between 0 and 1")


def _choice(name: str, value: object, allowed: tuple[str, ...]) -> None:
    if not isinstance(value, str) or value not in allowed:
        raise ValueError(f"{name} must be one of {allowed!r}")


# Multiple-candidate corrections for the smear significance test (hybrid.lengths):
# "bonferroni" multiplies each p value by the number of below-top candidates tested.
SMEAR_CORRECTIONS = ("bonferroni", "none")


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
        _number(
            "allele_selection.secondary_mode_min_fraction",
            self.secondary_mode_min_fraction,
            0.0,
            1.0,
        )
        if self.secondary_mode_min_fraction == 0:
            raise ValueError("allele_selection.secondary_mode_min_fraction must be > 0")
        _integer("allele_selection.min_allele_primary_records", self.min_allele_primary_records, 1)


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
class HybridSettings:
    """Read-centric engine thresholds (spec 2026-09-23 §5).

    Every default is provisional (prototype-derived) and is tuned on the benchmark
    dev/validation splits only; the sealed test split never informs a default.
    """

    anchor_max_edits: int = 12
    min_span_units: int = 15
    max_span_units: int = 160
    peak_window_base_bp: float = 30.0
    peak_window_per_unit_bp: float = 0.6
    min_peak_reads: int = 8
    far_peak_min_frac: float = 0.03
    near_peak_min_frac: float = 0.20
    rejected_peak_noise_reads: int = 2
    n_poa: int = 40
    poa_backend: str = "pyabpoa"
    polish_rounds: int = 2
    hp_vote: bool = True
    # S3/S7 (Task 6) POA sampling window and pileup/homopolymer-vote polishing tunables.
    poa_sample_window_floor_bp: float = 15.0
    poa_sample_window_frac: float = 0.006
    polish_insertion_majority_frac: float = 0.5
    hp_vote_min_run: int = 4
    het_af_min: float = 0.2
    het_min_group: float = 0.15
    link_phi_min: float = 0.5
    min_linked_sites: int = 2
    min_fragment_bp: int = 1000
    assign_margin: int = 3
    assign_max_error_rate: float = 0.15
    # S5/S6 (Task 7): ladder-flank width wrapped around each allele draft to build the
    # references that reads are assigned against.
    assign_flank_bp: int = 500
    qc_residual_af: float = 0.25
    max_unassigned_spanning_fraction: float = 0.2
    depth_adequate_spanning: int = 30
    depth_low_spanning: int = 10
    hp_llr_min: float = 10.0
    hp_min_reads: int = 20
    hp_min_alt_frac: float = 0.30
    hp_min_strand_reads: int = 5
    seed: int = 1
    # S1 anchor-search tunables (flank-anchor fallback when a motif is mutated).
    flank_anchor_bp: int = 30
    flank_anchor_edit_divisor: int = 4
    flank_anchor_edit_floor: int = 2
    # S2 length-model tunables: KDE shape, peak spacing, and the smear/rejection rules.
    kde_bandwidth_base_bp: float = 8.0
    kde_bandwidth_per_bp: float = 0.004
    kde_kernel_truncation_bw: float = 4.0
    kde_grid_step_bp: float = 2.0
    kde_grid_margin_bp: float = 100.0
    smear_short_product_units: float = 1.5
    peak_far_near_boundary_units: float = 2.0
    peak_min_separation_units: float = 0.7
    # C4.2 fix round 4: a below-top candidate is smear unless its read count in a core
    # window (smear_test_window_frac x its assignment half-window) significantly exceeds
    # the local smear background on BOTH sides (exact conditional Poisson rate test, one
    # sided, at smear_test_alpha after smear_test_correction over the candidates tested).
    # Adjusted p in [alpha / factor, alpha * factor) is the borderline band
    # (smear_ambiguous). Each background side starts at the assignment-window edge, spans
    # at least smear_background_flank_units repeat units and widens until it holds
    # smear_background_min_reads reads (or reaches the below-top region's edge).
    smear_test_alpha: float = 0.001
    smear_test_borderline_factor: float = 3.0
    smear_test_correction: str = "bonferroni"
    smear_test_window_frac: float = 0.25
    smear_background_flank_units: float = 2.0
    smear_background_min_reads: int = 5

    def __post_init__(self) -> None:
        for name in (
            "anchor_max_edits",
            "min_peak_reads",
            "rejected_peak_noise_reads",
            "polish_rounds",
            "min_linked_sites",
            "min_fragment_bp",
            "assign_margin",
            "depth_low_spanning",
            "hp_min_reads",
            "hp_min_strand_reads",
            "seed",
        ):
            _integer(f"hybrid.{name}", getattr(self, name))
        _integer("hybrid.n_poa", self.n_poa, 1)
        _integer("hybrid.hp_vote_min_run", self.hp_vote_min_run, 2)
        _integer("hybrid.min_span_units", self.min_span_units, 1)
        _integer("hybrid.max_span_units", self.max_span_units, self.min_span_units + 1)
        _integer(
            "hybrid.depth_adequate_spanning", self.depth_adequate_spanning, self.depth_low_spanning
        )
        _integer("hybrid.flank_anchor_bp", self.flank_anchor_bp, 1)
        _integer("hybrid.flank_anchor_edit_divisor", self.flank_anchor_edit_divisor, 1)
        _integer("hybrid.flank_anchor_edit_floor", self.flank_anchor_edit_floor, 0)
        _integer("hybrid.assign_flank_bp", self.assign_flank_bp, 1)
        _number("hybrid.peak_window_base_bp", self.peak_window_base_bp, 1)
        _number("hybrid.peak_window_per_unit_bp", self.peak_window_per_unit_bp)
        _number("hybrid.poa_sample_window_floor_bp", self.poa_sample_window_floor_bp, 0)
        _number("hybrid.poa_sample_window_frac", self.poa_sample_window_frac, 0)
        for name in (
            "far_peak_min_frac",
            "near_peak_min_frac",
            "het_min_group",
            "qc_residual_af",
            "hp_min_alt_frac",
            "link_phi_min",
            "assign_max_error_rate",
            "max_unassigned_spanning_fraction",
            "polish_insertion_majority_frac",
        ):
            _number(f"hybrid.{name}", getattr(self, name), 0, 1)
        _number("hybrid.het_af_min", self.het_af_min, 0.01, 0.5)
        _number("hybrid.hp_llr_min", self.hp_llr_min)
        _boolean("hybrid.hp_vote", self.hp_vote)
        _choice("hybrid.poa_backend", self.poa_backend, ("pyabpoa", "pyspoa"))
        _number("hybrid.kde_bandwidth_base_bp", self.kde_bandwidth_base_bp, 1.0)
        _number("hybrid.kde_bandwidth_per_bp", self.kde_bandwidth_per_bp, 0)
        _number("hybrid.kde_kernel_truncation_bw", self.kde_kernel_truncation_bw, 1.0)
        _number("hybrid.kde_grid_step_bp", self.kde_grid_step_bp, 0.1)
        _number("hybrid.kde_grid_margin_bp", self.kde_grid_margin_bp, 0)
        _number("hybrid.smear_short_product_units", self.smear_short_product_units, 0.01)
        _number("hybrid.peak_far_near_boundary_units", self.peak_far_near_boundary_units, 0)
        _number("hybrid.peak_min_separation_units", self.peak_min_separation_units, 0)
        _open_unit_interval("hybrid.smear_test_alpha", self.smear_test_alpha)
        _number("hybrid.smear_test_borderline_factor", self.smear_test_borderline_factor, 1)
        _choice("hybrid.smear_test_correction", self.smear_test_correction, SMEAR_CORRECTIONS)
        _number("hybrid.smear_test_window_frac", self.smear_test_window_frac, 0, 1)
        if self.smear_test_window_frac == 0:
            raise ValueError("hybrid.smear_test_window_frac must be > 0")
        _number("hybrid.smear_background_flank_units", self.smear_background_flank_units, 0)
        if self.smear_background_flank_units == 0:
            raise ValueError("hybrid.smear_background_flank_units must be > 0")
        _integer("hybrid.smear_background_min_reads", self.smear_background_min_reads, 1)


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
