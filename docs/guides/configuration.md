# Runtime configuration

MucOneSpan accepts a versioned JSON configuration for execution defaults,
repeat-selection and classification heuristics, consensus boundaries, confidence
weights, and reference layout. Defaults come from immutable typed settings in
`muc_one_span.settings`. The repository's `examples/runtime-settings.json` is a
complete default configuration generated with `muconespan settings show`.

## Use a configuration file

Put the global `--config` option before the command:

```bash
muconespan --config settings.json run \
  --input reads.bam \
  --output-dir results \
  --threads 8
```

Values are selected in this order:

1. An explicitly supplied command-line option.
2. The corresponding configuration-file value.
3. The central default.

Thus `--threads 8` overrides `run.threads` in the file. Omitted fields and sections
retain their central defaults. An explicit `--no-report` overrides
`"report": true`. Input and output locations remain command arguments; they are
not additional JSON fields. Use `muconespan COMMAND --help` for available flags.

A minimal configuration is:

```json
{
  "schema_version": 1,
  "run": {
    "threads": 4,
    "platform": "hifi"
  }
}
```

A supplied file must declare integer `schema_version: 1`. The root and each
section must be JSON objects. Unknown fields, duplicate keys, incorrect types,
nonfinite numbers, and invalid ranges are rejected. JSON booleans are `true` and
`false`; strings such as `"false"` are not accepted as booleans. Integer parameters
reject booleans and floating-point values such as `4.0`.

Relative `run.reference`, `run.clair3_model`, and `repeat_dictionary` paths resolve
against the configuration file's directory. Relative paths explicitly supplied
on the command line retain their normal command-line interpretation from the
working directory. `null` selects the bundled dictionary/reference or automatic
preset where supported. An empty `clair3_model` preserves the external caller's
model-selection behavior; it does not verify which model the caller uses.

For example, a file at `project/settings.json` containing
`"clair3_model": "models/hifi"` selects `project/models/hifi`. Select the model and
platform together. Executables continue to resolve through the existing tool
handling and `PATH`; configuration does not embed a particular installation path.

## Parameters

The defaults below describe the current configuration schema. They are preserved
scientific defaults, not recommendations to tune against a validation sample.

### Execution and allele selection

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `run.threads` | `4` | Per-tool thread setting; integer >=1, not a process-tree CPU cap. |
| `run.platform` | `"hifi"` | `hifi` or `ont`. |
| `run.min_coverage` | `10` | Minimum coverage used for allele selection; integer >=1. |
| `run.min_qual` | `5` | VCF QUAL filtering threshold; finite number >=0. |
| `run.clair3_model` | `""` | Model path; empty preserves the external default. |
| `run.reference` | `null` | Reference FASTA path; bundled ladder when omitted. |
| `run.minimap2_preset` | `null` | Explicit preset or automatic platform selection. |
| `run.mapping_timeout` | `3600.0` | Total mapping budget in seconds, including sorting, indexing and cleanup; finite number >0. `map`/`run --mapping-timeout` overrides it. |
| `run.report` | `false` | Generate the optional HTML report. |
| `allele_selection.min_gap` | `5` | Minimum gap between contig repeat counts separating clusters; integer >=1. |
| `allele_selection.valley_min_points` | `3` | Minimum points for valley consideration; integer >=3. |
| `allele_selection.valley_min_separation` | `3` | Minimum peak separation; integer >=1. |
| `allele_selection.refinement_max_shift` | `1` | Maximum refinement shift in repeat counts; integer >=0. |
| `allele_selection.secondary_mode_min_fraction` | `0.2` | Clinical gate. Allele selection is unresolved when a cluster's primary-alignment count at least `min_gap` units from its peak reaches this fraction of the peak. Number in (0,1]. Blocks a negative result; never creates a call. |
| `allele_selection.min_allele_primary_records` | `30` | Clinical gate. Minimum primary alignment records (a molecule proxy) per selected allele. Below it, `depth_status` is `low`, which blocks NEGATIVE and PATHOGENIC. Integer >=1. |
| `allele_selection.refinement_min_supported_records` | `3` | Peak refinement: minimum alignment records for a cluster contig to be a best-contig candidate. Formerly the literal `3` in `alleles.py` (`refine_peak_contig`). Integer >=1. |
| `allele_selection.refinement_supported_fraction` | `0.25` | Peak refinement: a candidate also needs this fraction of the best-covered contig's records (the larger of the two floors applies). Formerly the literal `0.25` in `alleles.py`. Number in (0,1]. |
| `allele_selection.refinement_min_shift_ont` | `2` | ONT only: the refined contig may move the cluster centre by up to the larger of this value and `refinement_max_shift`. Formerly the literal `2` in `alleles.py` (`_build_allele_info`). Integer >=0. |
| `allele_selection.valley_min_canonical_repeats` | `10` | Indel-valley splitter: valleys at or above this canonical repeat count are preferred when at least two exist. Formerly the literal `10` in `alleles.py` (`_split_cluster_by_indel`). Integer >=1. |
| `allele_selection.minority_min_alignment_records` | `3` | Minority-allele search: minimum alignment records per contig and per minority cluster. Formerly the literal `3` (both sites) in `alleles.py` (`detect_alleles`). Integer >=1. |
| `allele_selection.dominance_close_candidate_repeats` | `6` | Read dominance: a read aligned to only one of two candidates closer than this many repeats is ambiguous rather than dominant. Formerly the literal `6` in `read_dominance.py`. Integer >=1. |
| `allele_selection.dominance_zero_primary_extra_reads` | `2` | Read dominance: extra dominant reads required when the second candidate has no primary alignment. Formerly the literal `2` in `read_dominance.py`. Integer >=1. |
| `allele_selection.read_length_split_min_reads` | `5` | Read-length splitter: minimum reads per mode (twice this in total). Formerly the literal `5` in `length_candidates.py` (`split_cluster_by_read_length`). Integer >=1. |
| `allele_selection.read_length_split_min_fraction` | `0.15` | Read-length splitter: a mode needs this fraction of all reads. Formerly the literal `0.15` in `length_candidates.py`. Number in (0,1]. |
| `allele_selection.read_length_split_bin_bp` | `5` | Read-length splitter: histogram bin width in bases; a peak must dominate the two neighbouring bins on each side (structural). Formerly the literal `5` in `length_candidates.py`. Integer >=1. |
| `allele_selection.read_length_split_min_delta_bp` | `45` | Read-length splitter: minimum distance between the two modes in bases. Formerly the literal `45` in `length_candidates.py`. Integer >=1. |
| `allele_selection.read_length_split_max_delta_bp` | `320` | Read-length splitter: maximum distance between the two modes in bases. Formerly the literal `320` in `length_candidates.py`. Integer greater than `read_length_split_min_delta_bp`. |
| `allele_selection.read_length_split_unit_tolerance_bp` | `15` | Read-length splitter: the mode distance must lie within this many bases of a whole number of repeat units. Formerly the literals `15`/`45` in `length_candidates.py`. Integer >=0; at run time twice the value must be below the dictionary `repeat_length_bp`, otherwise a `ValueError` is raised. |
| `allele_selection.read_length_split_offset_bp` | `30` | Read-length splitter: non-repeat bases subtracted from a mode before conversion to repeat units. Formerly the literal `30` in `length_candidates.py`. Integer >=0. The unit length (formerly `60`) now comes from the repeat dictionary and the fixed repeat count (formerly `9`) from `reference_layout`. |

### Classification and consensus

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `classification.max_indel_probe` | `30` | Maximum indel-length search extent in bases; integer >=0. |
| `classification.max_fit_edit_distance` | `3` | Repeat-fitting edit-distance cutoff; integer >=0. |
| `classification.novel_repeat_edit_distance` | `2` | Distance above which a fitted unit is labeled novel; integer >=0. |
| `classification.minimum_unit_fraction` | `0.5` | Minimum fraction of nominal unit length considered; finite number in (0,1]. |
| `classification.early_stop_edit_distance` | `1` | Search early-stop threshold; integer >=0. |
| `classification.strict_segmentation` | `false` | Apply stricter fitting and unresolved-sequence handling. |
| `consensus.flank_length` | `500` | Reference flank extent per side, in bases; integer >=0. |
| `consensus.anchor_bases` | `20` | Up to this many available bases from each side of an exact boundary anchor; integer >=1. |
| `consensus.anchor_tolerance` | `50` | Inclusive displacement from the expected VNTR boundary in bases; integer >=0. |

`consensus.haploid_majority` and `consensus.haploid_min_qual` are deprecated and
have no effect. Configurations that set non-default values load with a logged
deprecation warning, which the CLI shows; use the `calling.*` fields.

### Confidence weights

Confidence values describe heuristic dictionary fit and VCF concordance. They
are not calibrated probabilities, and changing a weight cannot establish exact
variant support. Sequence-concordance requirements remain unchanged.

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `confidence.qual_low` | `5` | Lower QUAL breakpoint; finite number >=0. |
| `confidence.qual_high` | `20` | Upper QUAL breakpoint; finite number greater than `qual_low`. |
| `confidence.weight_below` | `0.3` | Weight below the lower breakpoint. |
| `confidence.weight_low` | `0.5` | Weight at the lower breakpoint. |
| `confidence.weight_high` | `1` | Weight at or above the upper breakpoint. |
| `confidence.absent_weight` | `0.3` | Weight for verified absence of matching VCF support. |
| `confidence.boundary_repeats` | `3` | Terminal-repeat region considered by the boundary rule; integer >=0. |
| `confidence.boundary_penalty` | `0.5` | Multiplier where the boundary penalty applies. |

All weights and the boundary multiplier must be finite numbers in [0,1]. Between
the QUAL breakpoints, the weight is linearly interpolated. `absent_weight` applies
only when projection succeeds and no matching VCF event is present. A
`projection_unavailable`, `localization_ambiguous` or
`heterozygous_genotype_unresolved` support status does not apply that absence
penalty and does not claim variant support. The independent boundary
rule can still apply `boundary_penalty` to a terminal repeat.

### Calling and optional read phasing

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `calling.sample_name` | `"sample"` | VCF sample name; nonempty, without whitespace or control characters. |
| `calling.read_phase` | `false` | Explicit opt-in to experimental same-length read phasing. |
| `calling.haploid_majority` | `true` | Apply the allele-fraction genotype rule to length-partitioned calls and to read-phased (haplotagged) haplotype calls. With `false`, diploid genotypes are kept on both paths, and any heterozygous record then leaves that allele unresolved (IUPAC candidate, no independent haplotype credit). |
| `calling.haploid_min_qual` | `4.0` | QUAL threshold for length-partitioned and read-phased haplotype calls; number >=0 or `null`. `null` applies `run.min_qual`/`--min-qual`. A non-default `--min-qual` that this value overrides is logged. The applied value is recorded in `alleles.json` as `variant_filter`. |
| `calling.haploid_alt_fraction` | `0.5` | Length-partitioned and read-phased haploid calls: the ALT fraction of allele-specific reads (`FORMAT/AD`) at or above which the genotype becomes ALT. Number in [0,1]. |
| `calling.haploid_ref_fraction` | `0.2` | Length-partitioned and read-phased haploid calls: ALT fraction below which the genotype becomes REF (`0/0`). Values in between keep the heterozygous call, which leaves that allele unresolved on both paths. Must be below `haploid_alt_fraction`. |
| `calling.stage_discordance_min_af` | `0.5` | Clair3 pileup FORMAT/AF (ALT reads / DP) at or above which a frameshift missing from the applied calls blocks NEGATIVE (`stage_concordance`). Chosen equal to `calling.haploid_alt_fraction`, the partition's haploid ALT cut-off. Measured: missed frameshifts at AF 0.68-0.85; the highest AF among unflagged NEGATIVE results was 0.156. Number in (0, 1]. |
| `calling.stage_discordance_min_depth` | `10` | Minimum pileup FORMAT/DP for such a record. Chosen equal to the `run.min_coverage` default. Measured: missed frameshifts at DP 33-408; a flank artefact at DP 2 must be excluded, and every value from 3 to 33 gave identical decisions. Integer >=1; a lower value is the conservative direction for this veto. |
| `read_phasing.internal_downsampling` | `null` | Optional WhatsHap internal downsampling override; integer >=1. |
| `read_phasing.mapping_quality` | `null` | Optional WhatsHap mapping-quality override; integer >=0. |
| `read_phasing.min_haplotype_reads` | `5` | Experimental read-phased path: minimum reads per haplotype for the haplotag split. Formerly `min_dp` (fixed at `5`) in `calling.py`. Integer >=1. |

Read phasing remains disabled by default. There is no dedicated read-phasing CLI
flag; an explicit configuration can enable `calling.read_phase`. Enabling it is
an experimental choice, not a claim of improved accuracy. It does not recover
reads excluded by upstream allele selection or prove that a collapsed cluster
contains independently recovered diploid haplotypes.

When an override is `null`, MucOneSpan leaves that WhatsHap parameter unspecified
and uses the installed tool's default. When the helper invokes WhatsHap, its
metadata records the exact argument vector under `phase_command`, along with the
version and overrides. Missing tooling and
unresolved phase evidence have explicit statuses; execution failures from an
available tool remain visible. The helper preserves its genotype and common
phase-set checks. It does not introduce a read-count support floor.

### Clinical decision

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `clinical_decision.max_ambiguous_bases` | `10` | Summed classification `ambiguous_bases` across both alleles above which the report banner becomes INCONCLUSIVE; integer >=0. |
| `clinical_decision.legacy_min_total_reads` | `30` | Total-read fallback for a summary without per-allele `depth_status` (see the per-allele gates above); below it, the banner becomes INCONCLUSIVE. Integer >=1. |

These thresholds are resolved by `report.compute_clinical_decision`: an
explicit `settings` argument overrides a `clinical_decision` section recorded
under a summary's `configuration.settings` (written by a prior `run`), which
overrides these central defaults. The resolved values and their source
(`explicit`, `recorded_configuration`, or `default`) are returned under the
decision's `"thresholds"` key. A recorded section with an unknown field or an
out-of-range value raises the same `ValueError` as an invalid configuration
file; it is never silently ignored or defaulted.

### Dictionary, selected layout and ladder range

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `repeat_dictionary` | `null` | Custom repeat dictionary JSON path or bundled resource. |
| `reference_layout.pre` | `["1","2","3","4","5"]` | Ordered fixed repeat IDs before the variable region. |
| `reference_layout.after` | `["6","7","8","9"]` | Ordered fixed repeat IDs after the variable region. |
| `reference_layout.min_units` | `1` | Minimum variable-region canonical repeat count in a generated ladder; integer >=1. |
| `reference_layout.max_units` | `150` | Maximum variable-region canonical repeat count; integer >=`min_units`. |

Each layout list must be nonempty, contain distinct nonempty strings, and be
disjoint from the other list. Selected IDs must exist in the loaded dictionary.
The fixed repeat count is `len(pre) + len(after)`; outer boundary anchors use the
first selected `pre` ID and last selected `after` ID. Dictionary category lists
can contain alternatives and do not define this selected layout.

The ladder range counts variable canonical repeats. A contig with N variable
repeats has N plus the selected fixed repeat count total VNTR repeats. The
`min_units` and `max_units` settings control ladder generation; they do not
regenerate or restrict an already supplied reference.

Generate a ladder using the same configuration intended for analysis:

```bash
muconespan --config settings.json ladder --output matching_ladder.fa
muconespan --config settings.json run \
  --input reads.bam --reference matching_ladder.fa --output-dir results
```

Changing the dictionary, selected layout or flank length requires a matching
explicit reference for `run`. Supply it through `run.reference` or `--reference`.
The application cannot infer that an arbitrary FASTA was built from matching
settings; keep the generating configuration with that reference.

## Hybrid Engine (Experimental)

`muconespan run --engine hybrid` selects a read-centric allele reconstruction
path (motif anchoring, a length model, partial-order-alignment consensus,
linked-site phase splitting, all-read assignment, polishing, and per-event
read-level support) instead of the default ladder-alignment/Clair3 path. See
[Core Concepts](../getting-started/concepts.md#hybrid-engine-experimental) for
the stage-by-stage pipeline and
[Known Limitations](../reference/limitations.md#hybrid-engine-experimental)
for measured detection limits and validation numbers.

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --engine hybrid \
  --assay amplicon \
  --threads 8
```

`--engine` is `ladder` or `hybrid` (`run.engine`); `--assay` is `amplicon` or
`genomic` (`run.assay`) and is recorded for provenance
(`summary["hybrid"]["assay"]`) -- it does not currently change any
`hybrid.*` default. `--report-igv` is rejected with `--engine hybrid` (a
hybrid run has no BAM alignment tracks to show). `--clair3-model` and
`--min-qual` are accepted but unused by the hybrid path.

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `run.engine` | `"ladder"` | `"ladder"` or `"hybrid"`. `ladder` stays the default until the benchmark decision rule is met on the sealed test split. |
| `run.assay` | `"amplicon"` | `"amplicon"` or `"genomic"`; library type, recorded only. |

### Install the `hybrid` extra

```bash
pip install 'muc_one_span[hybrid]'
```

installs `edlib`, `pyabpoa`, and `pyspoa`. Without the extra, `--engine
hybrid` fails with a clear `ImportError` naming the extra
(`"The hybrid engine needs the 'hybrid' extra: pip install
'muc_one_span[hybrid]'"`); it never silently falls back to another backend.
All three packages are MIT-licensed; the hybrid engine does not use medaka
or dorado.

| Package | Wheels | Notes |
| --- | --- | --- |
| `edlib` | manylinux/musllinux/macOS wheels on 3.10-3.13 | No wheel on 3.14 yet; the sdist builds and imports from source with a C compiler. |
| `pyabpoa` (default backend) | **sdist only** | Always builds from source; needs a C compiler and zlib (`gcc`, `libc6-dev`, `zlib1g-dev` on Debian/Ubuntu). Bioconda ships binaries. |
| `pyspoa` (alternative backend, `hybrid.poa_backend: "pyspoa"`) | manylinux wheels (x86_64, aarch64) | **No macOS wheel**; the sdist needs cmake and a C++ compiler. |

The project's own Docker image installs `gcc`, `libc6-dev` and `zlib1g-dev`
in the builder stage to build `pyabpoa`; only the built virtual environment
is copied into the runtime image.

### Hybrid settings (`hybrid.*`)

Every default below is **provisional** (prototype-derived) and tuned on the
development/validation splits only; the sealed test split never informs a
default. Every threshold is a validated `HybridSettings` field -- there are
no hardcoded thresholds in the hybrid engine.

#### Anchoring and span categorization (S1)

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.anchor_max_edits` | `12` | Maximum edlib edit distance for a motif-1/motif-9 anchor (both strands); integer >=0. |
| `hybrid.min_span_units` | `15` | Minimum accepted spanning-read length, in repeat units; integer >=1. |
| `hybrid.max_span_units` | `160` | Maximum accepted spanning-read length, in repeat units; integer >= `min_span_units` + 1. |
| `hybrid.flank_anchor_bp` | `30` | Ladder flank length used as a fallback anchor when a motif anchor cannot be found; integer >=1. |
| `hybrid.flank_anchor_edit_divisor` | `4` | Divides `anchor_max_edits` to derive the flank-anchor edit budget; integer >=1. |
| `hybrid.flank_anchor_edit_floor` | `2` | Minimum flank-anchor edit budget (`max(floor, anchor_max_edits // divisor)`); integer >=0. |

#### Length model (S2)

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.peak_window_base_bp` | `30.0` | Assignment half-window base width, in bp; number >=1. |
| `hybrid.peak_window_per_unit_bp` | `0.6` | Extra half-window width per repeat unit of length (models span-length noise growing with length); number >=0. |
| `hybrid.min_peak_reads` | `8` | Absolute minimum reads a candidate length peak needs; integer >=0. |
| `hybrid.far_peak_min_frac` | `0.03` | Minimum support fraction of total spanning reads for a peak >= `peak_far_near_boundary_units` from the top peak; number in [0,1]. |
| `hybrid.near_peak_min_frac` | `0.20` | Minimum support fraction for a peak nearer than `peak_far_near_boundary_units`; number in [0,1]. |
| `hybrid.rejected_peak_noise_reads` | `2` | A candidate peak with at most this many reads is `noise`, not gate-relevant; integer >=0. |
| `hybrid.kde_bandwidth_base_bp` | `8.0` | Gaussian KDE bandwidth base, in bp; number >=1.0. |
| `hybrid.kde_bandwidth_per_bp` | `0.004` | Extra KDE bandwidth per bp of length; number >=0. |
| `hybrid.kde_kernel_truncation_bw` | `4.0` | KDE kernel truncation, in bandwidths; number >=1.0. |
| `hybrid.kde_grid_step_bp` | `2.0` | KDE evaluation grid step, in bp; number >=0.1. |
| `hybrid.kde_grid_margin_bp` | `100.0` | KDE grid margin beyond the observed length range, in bp; number >=0. |
| `hybrid.smear_short_product_units` | `1.5` | Below this many units under the top peak, a read is a "short product", never tested as an allele candidate; number >=0.01. |
| `hybrid.peak_far_near_boundary_units` | `2.0` | Distance from the top peak, in units, beyond which a candidate uses the "far" support fraction; number >=0. |
| `hybrid.peak_min_separation_units` | `0.7` | Minimum separation, in units, between kept KDE maxima; number >=0. |
| `hybrid.smear_test_alpha` | `0.001` | Significance level of the one-sided exact conditional Poisson smear test; number strictly in (0,1). |
| `hybrid.smear_test_borderline_factor` | `3.0` | Width of the borderline p-value band around alpha (`[alpha/factor, alpha*factor)`), reported as `smear_ambiguous` instead of a silent `smear` rejection; number >=1. |
| `hybrid.smear_test_correction` | `"bonferroni"` | Multiple-testing correction across below-top candidates tested; `"bonferroni"` or `"none"`. |
| `hybrid.smear_test_window_frac` | `0.25` | Core window width, as a fraction of the candidate's assignment half-window, scored against the local background; number strictly >0, in [0,1]. |
| `hybrid.smear_background_flank_units` | `2.0` | Minimum span, in units, of each side of the smear-test background window; number strictly >0. |
| `hybrid.smear_background_min_reads` | `5` | Minimum reads required in each widened background window side; integer >=1. |

#### POA draft and polishing (S3, S7)

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.n_poa` | `40` | Maximum spanning reads sampled (with the seeded RNG) for the POA draft; integer >=1. |
| `hybrid.poa_backend` | `"pyabpoa"` | `"pyabpoa"` or `"pyspoa"`; no silent fallback when the selected backend is unavailable. |
| `hybrid.polish_rounds` | `2` | Pileup + homopolymer-vote polishing rounds; integer >=0. |
| `hybrid.hp_vote` | `true` | Run the homopolymer median-length vote after each pileup round; boolean. |
| `hybrid.poa_sample_window_floor_bp` | `15.0` | Floor of the "near-modal" length window POA draft members are sampled from, in bp; number >=0. |
| `hybrid.poa_sample_window_frac` | `0.006` | Fraction of the median length added to the POA sampling window (`max(floor, frac * median)`); number >=0. |
| `hybrid.polish_insertion_majority_frac` | `0.5` | An insertion slot is accepted only when its winning vote exceeds this fraction of covering votes; number in [0,1]. |
| `hybrid.hp_vote_min_run` | `4` | Minimum consensus run length rewritten by the homopolymer median vote; integer >=2. |

#### Phase split (S4)

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.het_af_min` | `0.2` | Minimum allele fraction for a candidate phase site (run-length sites use `max(het_af_min, phase_run_bg_multiplier * background)`); number in [0.01, 0.5]. A minor allele below this floor never forms a candidate: at the default, an **equal-length heterozygote with a minor allele at 15-20% produces a silent `none` split, not a flag** (`het_af_min` > `het_min_group`, so a group at the `het_min_group` edge can never form). |
| `hybrid.het_min_group` | `0.15` | Minimum fraction of members the smaller phase group must reach, else `unconfirmed_group_size`; number in [0,1]. |
| `hybrid.link_phi_min` | `0.5` | Minimum absolute phi correlation for two candidate sites to be linked; number in [0,1]. |
| `hybrid.min_linked_sites` | `2` | Minimum linked events required to split a peak; fewer produces `unconfirmed_single_site`; integer >=0. |
| `hybrid.phase_max_site_reads` | `300` | Read cap for building the phase site table (sampled with the seeded RNG above the cap); integer >=1. |
| `hybrid.phase_run_min_len` | `3` | Minimum homopolymer run length treated as a run-length candidate site; integer >=2. |
| `hybrid.phase_run_bg_window` | `3` | Run-length background window (+/- d observed length) used to score a run site; integer >=1. |
| `hybrid.phase_min_minor_reads` | `5` | Minimum reads showing a candidate site's minor allele; integer >=1. |
| `hybrid.phase_run_bg_multiplier` | `4.0` | Multiplier on local background noise for the run-length candidate-site AF floor; number >=0. |
| `hybrid.phase_run_safety_multiplier` | `2.0` | NEGATIVE-blocking safety tier for an equal-length genotype (Task 15g). When the length model finds a **single** peak and that peak stays unsplit with no candidate site (`none`), every homopolymer run is tested against the lower floor `max(het_af_min, phase_run_safety_multiplier * background)`, with the same leave-one-out peer background as a candidate site but none of its other tests. A run above it gives the phase basis `unconfirmed_run_site` (selection and phase status `unresolved_run_site`): the result is INCONCLUSIVE with the located reason `unresolved heterozygous site at repeat N`. The tier never splits the peak, never creates an event and never makes a result PATHOGENIC. It exists because a heterozygous run can sit below the split floor: simulated HiFi reads of an equal-length heterozygous dupC showed 42% C8 against 13.7% at the peer C7 runs (ratio 3.1, under the x4 split floor), and were reported NEGATIVE before this tier. Default from the v4 development panels (every wild-type length peak left unsplit, whatever the peak count): the floor is reached at 11/26 HiFi peaks at 1.5, 4/26 at 2.0 and 3/26 at 2.5-3.0, and at none of 61 ONT peaks; 2.0 is the lowest value on that plateau. The wild-type runs that reach it carry 21-30% of reads at one length, at or above `het_af_min`; the tier itself flags 3 of the 4 equal-length wild-type HiFi samples and none of the 3 ONT ones. Setting it to `phase_run_bg_multiplier` restores the pre-15g behaviour for runs below the split floor. Number >=1. |
| `hybrid.phase_gap_af_factor` | `1.5` | AF factor applied when the candidate site's minor allele is a gap (deletion); number >=1. |
| `hybrid.phase_min_pair_reads` | `10` | Minimum reads informative at both sites of a pair before their linkage is tested; integer >=2. |
| `hybrid.phase_strand_bias_alpha` | `0.001` | Strand-bias test significance level; a site failing it, or with no minor-allele observation on a strand with >= `hp_min_strand_reads` reads, is rejected. Column and insertion sites use a one-sided Fisher exact test. Homopolymer-run sites use a stutter-aware test: each strand's length-error profile comes from the other runs of the same base, the minor run length's stutter-deconvolved share must reach `het_af_min`, and a likelihood-ratio test (one weight for both strands vs. one per strand) is applied at this level, so strand-asymmetric ONT stutter is not read as strand bias. Number strictly in (0,1). |
| `hybrid.phase_single_event_split` | `"indel"` | When the length model finds a **single** peak (an equal-length genotype) and that peak's candidate sites form exactly one heterozygous event (fewer linked events than `min_linked_sites`), split the peak on that event: `"indel"` only when the event changes the sequence length (every frameshift), `"all"` for any event, `"off"` never. Reads are grouped by their allele at the event (run sites by the more likely run length under the strand's stutter profile); the split is refused when either group is below `het_min_group` or both group drafts are identical. An unsplit peak stays `unconfirmed_single_site`, which blocks a negative call and names the site (`unresolved heterozygous site at repeat N` in `selection_detail`). With two length peaks a single event never splits a peak. The split also needs the peak-level share gate (`phase_single_event_alpha`). One of `"off"`, `"indel"`, `"all"`. |
| `hybrid.phase_single_event_alpha` | `0.001` | A single-event split selects each allele's reads by the event, so the event's read support is conditional on the split. The split is made only when the one-sided lower confidence bound (at this level; profile likelihood) of the stutter-deconvolved minor share reaches `het_af_min`, computed over a fresh sample of at most `phase_single_event_bound_reads` reads (drawn with `seed`, independent of the site table). The fixed sample keeps the bound's power independent of depth, since a systematic artefact does not shrink with depth. Otherwise the peak keeps the phase basis `unconfirmed_single_site` (selection status `unresolved_single_site`: INCONCLUSIVE, located). Number strictly in (0,1). |
| `hybrid.phase_single_event_bound_reads` | `300` | Sample size of the single-event share bound (`phase_single_event_alpha`). It is separate from the `phase_max_site_reads` compute cap, so raising that cap cannot make the bound more permissive; a larger value gives the bound more power at high depth. The default equals the former shared value. Integer >=1. |
| `hybrid.phase_run_error_cap` | `16` | Run-length error profiles of the stutter-aware run-site test pool errors beyond +/- this many bases into their edge bins; 16 (= `hp_max_run_len`) keeps every modelled run's full error range distinct. Integer >=1. |

#### Reference and read assignment (S5, S6)

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.assign_flank_bp` | `500` | Ladder flank width wrapped around each allele draft to build the reference every read is assigned against; integer >=1. |
| `hybrid.assign_margin` | `3` | Minimum edit-distance gap to the second-best reference before a read is assigned (else `undecided`); integer >=0. |
| `hybrid.assign_max_error_rate` | `0.15` | Reads needing more than this fraction of edits even to the best reference are `off_target`; number in [0,1]. |
| `hybrid.min_fragment_bp` | `1000` | Minimum length of a left/right-anchored or internal fragment considered for assignment; integer >=0. |

#### Depth thresholds

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.depth_adequate_spanning` | `30` | Spanning-read count at/above which `depth_status` is `adequate`; integer >= `depth_low_spanning`. |
| `hybrid.depth_low_spanning` | `10` | Spanning-read count at/above which `depth_status` is `low` (below it, `insufficient`); integer >=0. |

#### Residual QC and event read support (S8, S10)

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.qc_residual_af` | `0.25` | Minor-allele fraction at a consensus column (outside long runs) that becomes a `residual_heterogeneity` site; number in [0,1]. |
| `hybrid.qc_residual_min_run` | `3` | Consensus runs at/above this length are skipped by residual QC (a run indel has no unique column); integer >=2. |
| `hybrid.hp_event_min_run` | `4` | Minimum consensus run length for a single-base-indel dictionary template to be typed a homopolymer event (else it falls back to parent-vs-template competition); integer >=2. |
| `hybrid.hp_max_run_len` | `16` | Runs whose event or reference length would reach this cap are not modelled as a homopolymer mixture; integer >= `hp_event_min_run` + 1. |
| `hybrid.hp_background_pseudocount` | `0.5` | Additive smoothing pseudocount for the per-strand background run-length profile; number strictly >0. |
| `hybrid.hp_stutter_model` | `"length"` | Stutter background of the homopolymer event/no-event mixture. `"length"`: each allele is convolved with the per-strand stutter profile of **its own** run length (a dupC C8 run with the C8 profile, the no-event C7 allele with the C7 profile), measured at the sample's own peer runs of that base and length (the event run left out); a length without enough peer runs is extrapolated from the two nearest measured lengths of the base (no pooling across bases); with one measured length its profile is shifted; with none, the `"shift"` rule applies. `"shift"`: the no-event length's profile, shifted by the event for the event allele (behaviour before this setting). One of `"length"`, `"shift"`. |
| `hybrid.hp_stutter_min_class_runs` | `2` | Peer runs (same base and length, event run excluded) a length needs to be measured rather than extrapolated. One run is one sequence context, which cannot separate the effect of length from that of context. Integer >=1. |
| `hybrid.hp_stutter_min_class_reads` | `200` | Clean run observations on a strand a length needs to be measured on that strand. A length with enough peer runs but fewer observations is treated as unmeasured on that strand (extrapolated, guarded), because a profile from so few reads is dominated by the pseudocount and would make the alleles harder to tell apart. Integer >=1. |
| `hybrid.hp_stutter_max_growth` | `3.0` | Cap on the per-base growth (and, inverted, the shrinkage) of each error value's share when a length is extrapolated log-linearly from the two nearest measured lengths. The default sits above the largest growth between adjacent measured lengths in the development panels (about 2.7, C6 to C7 deletion). `1.0` disables growth (the nearest profile shifted). Number >=1. |
| `hybrid.hp_stutter_max_event_confusion` | `0.3` | Identifiability guard: an event-allele profile that is not measured (extrapolated or moved from another length) may put at most this share of its mass on the no-event run length; above it, that strand falls back to the `"shift"` rule. Where stutter saturates (real ONT "+" reads: deletion 10% at C6, 26% at C7 but 21% at a true C8) extrapolation would predict a C8 allele read as C7 about as often as C8, and a wild-type/dupC mixture could pass as a pure dupC. Development data: 0.03-0.19 on identifiable strands (simulated HiFi 0.17-0.19), 0.42-0.55 on saturating ONT "+" strands; the default lies between. `1.0` disables the guard. Number in [0,1]. |
| `hybrid.hp_llr_min` | `10.0` | Minimum stutter-aware log-likelihood ratio for a homopolymer event; number strictly >0. |
| `hybrid.hp_min_reads` | `20` | Minimum reads (`n`) before a status other than `insufficient_depth` is possible; integer >=1. |
| `hybrid.hp_min_alt_frac` | `0.30` | Minimum alt-supporting fraction for a homopolymer or competition event; number in [0,1]. |
| `hybrid.hp_min_strand_reads` | `5` | A strand with at least this many reads must not show a negative homopolymer LLR, else the event is `discordant`; integer >=0. |
| `hybrid.event_context_units` | `1.0` | Context, in repeat units (x the dictionary unit length), compared around a competition event's unit on each side; number >=0. |
| `hybrid.event_max_alternative_frac` | `0.25` | Maximum estimated alternative share at the event site: `ref/n` for competition events, `1 - f_hat` of the event/no-event stutter mixture (length-aware, `hp_stutter_model`) for homopolymer events; above it the event is `discordant`; number in [0,1]. |

#### Engine orchestration

| Section.field | Default | Meaning and validation |
| --- | --- | --- |
| `hybrid.polish_max_reads` | `120` | Maximum spanning+partial members sampled per allele for polishing; integer >=1. |
| `hybrid.polish_partial_min_units` | `1.0` | Minimum trimmed length, in repeat units, for an assigned non-spanning fragment to join the polishing pileup; number >=0. |
| `hybrid.qc_residual_max_reads` | `200` | Maximum spanning members sampled per allele for residual QC; integer >=1. |
| `hybrid.max_unassigned_spanning_fraction` | `0.2` | Above this fraction of spanning reads assigned to no allele, sample `selection_status` becomes `unresolved_unassigned_spanning`; number in [0,1]. |
| `hybrid.seed` | `1` | Seed for `random.Random` used by every random choice in the engine (POA/phase-table/reassignment sampling); deterministic given the same reads and settings; integer >=0. |

### Evidence fields

Per-allele, `depth_status` is `adequate`/`low`/`insufficient` from
`spanning_reads` alone (the genomic assay currently uses the **same**
thresholds as amplicon; `--assay genomic` records the library type but does
not lower them, so a low-molecule WGS run can legitimately show `low` or
`insufficient` depth even when the pipeline behaves correctly). Sample
`selection_status` is `"resolved"` or one of six `"unresolved_*"` reasons
(`unresolved_max_alleles`, `unresolved_single_site`,
**`unresolved_group_size`** -- a linked-site split whose smaller group falls
below `het_min_group` -- **`unresolved_run_site`** -- an unsplit single length
peak with a homopolymer run above the `phase_run_safety_multiplier` floor --
`unresolved_rejected_peak`, `unresolved_unassigned_spanning`); any `unresolved_*` selection blocks a
NEGATIVE result (`clinical_gates.allele_gate_reasons`) but does **not** block
PATHOGENIC when the causative event has its own explicit read-level support --
`compute_clinical_decision` only requires an unblocked mutation to reach
PATHOGENIC, and adds unresolved-selection reasons to that banner as "Quality
caveat" detail lines rather than withholding the call.

`"not_assessed"` is a **ladder-engine** `selection_status`/`depth_status`
value (`selection_qc.assess_allele`), used when alignment `fit_metrics` do
not carry enough information to judge selection or depth; the shared
`clinical_gates.depth_gate_failure` gate defers to the legacy total-read
fallback only while *no* allele in the sample carries an assessed depth
status. The hybrid engine always computes a per-allele `depth_status` from
`spanning_reads`, so it never emits `"not_assessed"` and that legacy
fallback never applies to a hybrid summary.

**Read-support evidence contract** (`hybrid/evidence.py`): only the spanning
reads assigned to the carrying allele, already oriented to its consensus,
count -- an unassigned or off-target read contributes nothing. An event is
`supported` only when those reads favour the event allele over *both* the
no-event allele (the unit reverted to its dictionary parent) and the
best read-derived alternative (the pileup/homopolymer-vote consensus of the
reads that do not favour the event), not just over one of the two; a tie
counts as neither. Homopolymer-run events (a single-base indel dictionary
template inside a consensus run >= `hp_event_min_run`) are fit with a
stutter-aware mixture instead of a raw vote: `event_allele_fraction` finds
the maximum-likelihood weight of the event allele in
`f * P(observed | event) + (1 - f) * P(observed | no-event)` over each
read's observed run length, using a per-strand background stutter profile
measured from the sample's other same-base, same-length runs; `1 - f_hat` is
`alternative_frac`, gated by `event_max_alternative_frac`. Status is one of
`supported`, `insufficient_depth` (`n < hp_min_reads`), `discordant`
(alternative share too high for either kind; for a homopolymer event, also
when a strand with >= `hp_min_strand_reads` reads shows a negative LLR),
`not_supported` (LLR or alt fraction below threshold), or `not_localized`
(the mutation's repeat unit or dictionary parent could not be found).

`consensus_concordance_fraction` (per allele) is the hybrid engine's own
read-support evidence: the mean, over every consensus position, of the
fraction of covering reads whose base agrees with the consensus, from the
same full/partial reads that built and polished that consensus
(`hybrid.polish.consensus_concordance`). It is reported alongside
`classification_confidence_status: "not_applicable_dictionary_fit_heuristic"`,
because the ladder's `classify.py` `confidence`/`allele_confidence` (the
dictionary-fit heuristic shown in `repeats.json`, the CLI's `confidence:`
line and the HTML report's "Allele confidence" tile) is computed identically
for both engines and carries no hybrid reconstruction evidence.

## Inspect and validate settings

`muconespan settings show` prints the effective settings as schema-1 JSON. The
output loads back with `--config`, so it is a complete starting point for a
custom file:

```bash
muconespan settings show > settings.json            # central defaults
muconespan settings show --config my.json           # defaults merged with my.json
muconespan settings show --section hybrid           # one section only
```

Without `--config` (either the command's own option or the global one placed
before `settings`), it prints the central defaults. `--config` runs the same
strict loader as `run`, so relative resource paths print as absolute paths
resolved against the file's directory. `--section NAME` prints only
`schema_version` and that section; the sections it omits keep their defaults
when the output is loaded. `examples/runtime-settings.json` is the output of
`muconespan settings show`, and a unit test keeps the two identical.

`muconespan settings validate FILE` runs the strict loader on `FILE` without
running anything else. A valid file prints a confirmation and exits 0. An
invalid file (unknown or duplicate keys, a missing `schema_version`, wrong
types, out-of-range values or malformed JSON) prints the loader's first error,
naming the file, and exits non-zero:

```bash
muconespan settings validate settings.json
```

To choose values for these settings from benchmark data rather than by hand,
see [Calibration](../benchmark.md#calibration). `benchsim calibrate` validates
every grid point with this same loader, and `calibrate-report` writes a
`recommended-config.json` that `--config` loads.

## Effective configuration and provenance

A `run` attempt records `run_configuration.json` in its output directory before
external tools execute, alongside the run-status record. It captures effective
settings after CLI overrides, the input configuration path/hash when supplied,
input, resolved-reference and repeat-dictionary paths/hashes, and model/tool-selection provenance.
Keep this file with `run_status.json`, tool versions, and analysis outputs when
comparing runs. Early configuration or file-validation failures can prevent a
complete input/reference hash record. A new invocation invalidates prior configuration
provenance before input checks; an output file alone does not establish a successful run.

Automatic preset selection remains null in the reusable settings block;
`resolved_minimap2_preset` records the actual preset used. Explicit resource/model
paths are made absolute in recorded settings. Changing platform on a later run
therefore retains automatic preset selection unless a preset was explicitly set.

Global configuration parsing occurs before the run callback. A malformed global
`--config` exits2 without touching any previous output directory; those old files
are not evidence of success for the failed invocation. Evaluation drivers must
record the command exit status, including failures before the callback begins.

This provenance describes the settings actually selected. It does not substitute
for recording the model files and installed-tool versions needed to reproduce a
run, particularly when external defaults were left implicit.

## Python use

Library callers can load the same file or construct validated immutable settings:

```python
from dataclasses import replace
from pathlib import Path

from muc_one_span.settings import DEFAULT_SETTINGS, RunSettings, RuntimeSettings, load_settings

settings = load_settings(Path("settings.json"))
operational = RuntimeSettings(run=RunSettings(threads=8))
adjusted = replace(DEFAULT_SETTINGS, run=replace(DEFAULT_SETTINGS.run, threads=8))
```

Pass the relevant settings section to a scientific function's `settings` keyword,
and `reference_layout` where accepted. Existing explicit function arguments take
precedence over their corresponding settings. Constructing or loading settings
does not mutate process-global state or invoke tools. Use `settings_as_dict` to
serialize an effective configuration; use `dataclasses.replace` to make a changed
copy rather than assigning to a frozen field.

Configuration covers operational choices and heuristic tunables. Codon size,
coordinate conversions, VCF genotype indexing, exact unit-cost edit-distance
recurrence and ties, supported-event identity, failure handling, and one-to-one
evaluation are algorithm or format contracts. They remain fixed. Evaluation
endpoints and acceptance tolerances also remain fixed by the validation protocol.
Nondefault settings require their own validation; this interface does not claim
scientific equivalence or calibrated confidence for arbitrary parameter choices.


The thread setting retains existing execution semantics. Mapping may run minimap2
and samtools sort concurrently, and two allele jobs can run concurrently; tool
thread flags can count additional workers differently. It is not a strict global
CPU budget. Compare measured CPU and wall times with the recorded configuration.

Standalone `consensus` now loads the bundled dictionary or configured/explicit
`--repeats-db` and uses the same anchor policy as `run`. This corrects its previous
fixed-only trimming. Anchors use the actual ladder-selected flank prefixes; their
components can be shorter than `anchor_bases`, and coordinates use actual lengths.
Zero tolerance accepts an anchor exactly at the expected boundary.

A run that enters the pipeline callback creates its output directory and execution
status before checking input existence. A misspelled input path therefore leaves
a failure record in that directory.
