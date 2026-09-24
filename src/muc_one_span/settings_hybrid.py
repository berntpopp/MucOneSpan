"""Read-centric hybrid-engine settings, split out of ``muc_one_span.settings``.

Kept separate so the hybrid engine's large, fast-growing tunable set does not
push ``settings.py`` past the repository's file-size gate. Imports only the
shared validation helpers from ``muc_one_span.settings_validation``; nothing
here imports back from ``muc_one_span.settings``.
"""

from __future__ import annotations

from dataclasses import dataclass

from muc_one_span.settings_validation import (
    _boolean,
    _choice,
    _integer,
    _number,
    _open_unit_interval,
)

# Multiple-candidate corrections for the smear significance test (hybrid.lengths):
# "bonferroni" multiplies each p value by the number of below-top candidates tested.
SMEAR_CORRECTIONS = ("bonferroni", "none")


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
    # S8/S10 (Task 9) read-level event evidence. A dictionary template that is a
    # single-base indel inside a consensus run >= hp_event_min_run is typed a
    # homopolymer event; runs are modelled up to hp_max_run_len (an event whose run
    # would reach that cap falls back to parent-vs-template competition), and the
    # per-strand background length profile adds hp_background_pseudocount per length.
    # Residual QC skips consensus runs >= qc_residual_min_run (column-ambiguous indels).
    hp_event_min_run: int = 4
    hp_max_run_len: int = 16
    hp_background_pseudocount: float = 0.5
    qc_residual_min_run: int = 3
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
    # S4 (Task 8) linked-site phase split: site-table read cap, homopolymer-run site
    # length and background window, minor-allele read floor, run background multiplier
    # (D6), gap-allele AF factor and pairwise-linkage read floor. Strand consistency is
    # a strand-bias test: a site is rejected when a one-sided Fisher exact test finds
    # its minor allele depleted on either strand at phase_strand_bias_alpha, or when
    # the minor is absent from a strand with >= hp_min_strand_reads reads.
    # phase_run_min_len (3) < hp_vote_min_run (4) on purpose (prototype values): the
    # split must treat a 3-run as one run-length site, since a run indel is ambiguous
    # per column; the polish median vote only needs to rewrite runs >= 4, where the
    # column-wise pileup vote is unreliable, and leaves shorter runs to that vote.
    phase_max_site_reads: int = 300
    phase_run_min_len: int = 3
    phase_run_bg_window: int = 3
    phase_min_minor_reads: int = 5
    phase_run_bg_multiplier: float = 4.0
    phase_gap_af_factor: float = 1.5
    phase_min_pair_reads: int = 10
    phase_strand_bias_alpha: float = 0.001
    # Engine orchestration (Task 11). Polishing and residual QC use at most
    # polish_max_reads / qc_residual_max_reads spanning members per allele (sampled with
    # the seeded RNG when a group is larger); an assigned non-spanning fragment joins the
    # polishing pileup only when its trimmed length reaches polish_partial_min_units
    # repeat units (unit length from the repeat dictionary).
    polish_max_reads: int = 120
    polish_partial_min_units: float = 1.0
    qc_residual_max_reads: int = 200

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
            "hp_min_strand_reads",
            "seed",
        ):
            _integer(f"hybrid.{name}", getattr(self, name))
        _integer("hybrid.n_poa", self.n_poa, 1)
        _integer("hybrid.hp_vote_min_run", self.hp_vote_min_run, 2)
        _integer("hybrid.hp_event_min_run", self.hp_event_min_run, 2)
        _integer("hybrid.hp_max_run_len", self.hp_max_run_len, self.hp_event_min_run + 1)
        _number("hybrid.hp_background_pseudocount", self.hp_background_pseudocount, 0)
        if self.hp_background_pseudocount == 0:
            raise ValueError("hybrid.hp_background_pseudocount must be > 0")
        _integer("hybrid.qc_residual_min_run", self.qc_residual_min_run, 2)
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
        _integer("hybrid.hp_min_reads", self.hp_min_reads, 1)
        _number("hybrid.hp_llr_min", self.hp_llr_min)
        if self.hp_llr_min == 0:
            raise ValueError("hybrid.hp_llr_min must be > 0")
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
        for name, minimum in (
            ("phase_max_site_reads", 1),
            ("phase_run_min_len", 2),
            ("phase_run_bg_window", 1),
            ("phase_min_minor_reads", 1),
            ("phase_min_pair_reads", 2),
        ):
            _integer(f"hybrid.{name}", getattr(self, name), minimum)
        _number("hybrid.phase_run_bg_multiplier", self.phase_run_bg_multiplier, 0)
        _number("hybrid.phase_gap_af_factor", self.phase_gap_af_factor, 1)
        _open_unit_interval("hybrid.phase_strand_bias_alpha", self.phase_strand_bias_alpha)
        _integer("hybrid.polish_max_reads", self.polish_max_reads, 1)
        _integer("hybrid.qc_residual_max_reads", self.qc_residual_max_reads, 1)
        _number("hybrid.polish_partial_min_units", self.polish_partial_min_units, 0)
