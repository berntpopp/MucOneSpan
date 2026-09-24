# Runtime configuration

MucOneSpan accepts a versioned JSON configuration for execution defaults,
repeat-selection and classification heuristics, consensus boundaries, confidence
weights, and reference layout. Defaults come from immutable typed settings in
`muc_one_span.settings`. The repository's `examples/runtime-settings.json` is a
complete default configuration generated from that API.

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
| `calling.stage_discordance_min_af` | `0.5` | Clair3 pileup FORMAT/AF (ALT reads / DP) at or above which a frameshift missing from the applied calls blocks NEGATIVE (`stage_concordance`). Equals `calling.haploid_alt_fraction`, the partition's haploid ALT cut-off. Measured: missed frameshifts at AF 0.68-0.85; the highest AF among unflagged NEGATIVE results was 0.156. Number in (0, 1]. |
| `calling.stage_discordance_min_depth` | `10` | Minimum pileup FORMAT/DP for such a record. Equals the `run.min_coverage` default. Measured: missed frameshifts at DP 33-408; a flank artefact at DP 2 must be excluded, and every value from 3 to 33 gave identical decisions. Integer >=1; a lower value is the conservative direction for this veto. |
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
