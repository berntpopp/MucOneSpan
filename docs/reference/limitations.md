# Known Limitations

## Allele Length Detection

### Discrete Repeat-Unit Length Errors

When length detection fails on called alleles, the error is almost exclusively in discrete integer multiples of 1 VNTR unit ($\pm 1$ repeat = $\pm 60$ bp, rarely $\pm 2$ to $3$ repeats). Because the reference ladder contigs are generated in discrete 60 bp steps, alignment peak finding and contig assignment naturally snap to repeat unit boundaries.

### Homozygous Same-Length Alleles and Read-Dominance Testing

When both alleles share the same repeat count (e.g. 60/60) or differ by a narrow gap (1–3 repeats), naive peak splitting or valley thresholding can over-split reads. In v0.12.0, pairwise read-dominance scoring (`read_dominance.py`) requires a calibrated affine-gap score margin ($\delta = 43$ HiFi, $42$ ONT) and a minimum count of dominant spanning reads ($\ge 3$ HiFi, $\ge 4$ ONT) before accepting a second allele. If the margin or support floor is not met, the pipeline conservatively reports a single candidate with multiplicity unresolved rather than forcing a false split.

### Extreme PCR Bias (Asymmetric Alleles and Amplicon Dropout)

In PCR amplicon mode with highly asymmetric allele pairs (e.g., 25/140, 20/90, 30/120), PCR exponential amplification heavily favors the shorter fragment. At typical template depths (e.g. 60 templates), the longer allele frequently yields only 0–1 usable reads. Because the pipeline enforces strict evidence floors to protect normal homozygous controls from false splitting, it correctly withholds calling the longer allele. In benchmark evaluations, an uncalled second allele appears as a missing prediction, producing an apparent error equal to the full allele gap (e.g. $\Delta = 115$ repeats / $6,900$ bp).

### Non-Spanning Read Limitations in Whole-Genome ONT Sequencing

In whole-genome sequencing (without targeted PCR amplicons), sequencing reads are randomly sheared across the genome with read-length distributions typically capped below 8–10 kb (e.g., 7.5 kb in standard NanoSim simulation profiles). Because real-world *MUC1* VNTR alleles span 3,000–7,000 bp and require several hundred base pairs of unique flank sequence on both sides to anchor assembly, the vast majority of whole-genome ONT reads start or terminate inside the VNTR array. Consequently, while whole-genome ONT demonstrates 100% normal control specificity (zero false alarms), its sensitivity to detect pathogenic variants in long alleles is constrained (~12–23% event recall) unless ultra-long reads (>20 kb) or targeted amplicon enrichment are employed.

On non-spanning whole-genome ONT reads, allele detection collapses reads across
the whole reference ladder into one cluster in 23/30 and 26/30 simulated
designs, and the reported length is more than 20 repeat units too long in
25/30 and 29/30 designs (0 false positives, 0 false negatives; every affected
case is reported INCONCLUSIVE, never a false PATHOGENIC or NEGATIVE). The
ladder caller is built and validated as a targeted amplicon caller. On
non-spanning genomic reads, its INCONCLUSIVE result on these designs is the
expected, correct outcome, not a defect: interpreting a genomic MUC1 VNTR
result requires either an amplicon assay or a purpose-built genomic mode,
neither of which this release provides.

## Variant Calling

### Long VNTR Alleles (>100 repeats)

Clair3 has reduced sensitivity for detecting 1bp insertions within very long tandem repeat regions (>6kb). For alleles with >100 repeat units, single-nucleotide frameshifts (e.g., dupC) may not be detected by Clair3, leading to false negatives.

**Affected test samples:**
- `sample_dupc_100_120`: dupC mutation in 100-repeat allele not detected
- `sample_long_120_140`: mutations in 120+ repeat alleles not detected

### Haploid allele-fraction rule without a depth floor

The AD-fraction haploid rule (`calling.haploid_alt_fraction`/`haploid_ref_fraction`) has no minimum informative depth yet, so a genotype on very few allele-specific reads can still be set to ALT or REF (tracked in #70).

### Boundary Repeat Mutations

Mutations detected in the last 3 repeat units of an allele receive a boundary penalty (0.5x confidence multiplier) because these positions are prone to alignment artifacts at the VNTR boundary.

## Classification

### Novel Repeat Types

Repeat sequences with >2 base substitutions from any known type are classified as `novel_repeat`. These may represent true biological variation or alignment/sequencing errors. The edit distance and identity percentage are reported for manual review.

### Unresolved segmentation

Poor-fitting windows and trailing sequence are reported explicitly. After an
uncertain window, repeat localization can be ambiguous even when the classifier
consumes the sequence. The previous gap-as-repeat recovery is not a reliable
biological reconstruction. An opt-in strict segmentation experiment retains
exact suffix evidence without inventing a repeat index, but is disabled by default
because it lost known events in development replay.

### Sequence vs. Length Discordance in Consensus Reconstruction

In diploid benchmarking, called alleles frequently achieve 100% exact repeat count (`count_exact: true`, $\Delta = 0$ bp) while failing strict nucleotide sequence identity (`sequence_exact: false`). The sequence edit distance is typically small (1 to 3 base pairs across the entire 3,000–6,000 bp array).

Empirical root-cause analysis across the simulation and validation cohorts reveals:

1. **IUPAC Ambiguity Codes from Clair3 `0/1` Calls (95.3% of mismatches):**
   The reference ladder contigs are composed of canonical `X` repeat units. Real alleles contain variant repeat units (`A`, `B`, `C`, `D`, etc.) that differ from `X` by 1–2 SNPs. When reads partitioned for an allele are aligned to the contig, repetitive cross-talk or sequencing errors lead to non-100% allele frequencies (e.g. 75–85% ALT). Clair3 frequently classifies these sites as heterozygous `0/1` rather than homozygous `1/1`. `bcftools consensus` then inserts IUPAC ambiguity codes (`S` for C/G, `M` for A/C, `R` for A/G, `Y` for C/T). While the overall length is 100% exact, IUPAC symbols count as mismatches in strict ACGT sequence comparisons, and prevent `classify.py` from matching any known pure ACGT repeat in the dictionary (rendering them as `?`).
2. **Catalogue Completeness (Zero Missing Units):**
   The repeat dictionary contains all 34 known biological units (fixed 1–9, canonical X, and variants A through W). Mismatches are not caused by missing units in the catalogue or random background SNPs from simulator tools.
3. **Clair3 False Negatives (Reference Fill):**
   In low-coverage regions, Clair3 occasionally misses a variant SNP entirely, leaving that repeat unit with the reference ladder's canonical `X` sequence instead of the true variant unit, resulting in a 1–2 bp substitution.
4. **Homopolymer Indels (ONT):**
   Nanopore reads occasionally suffer 1-bp indel compression within the 7-C homopolymer tract of unit `X`.


## ONT-Specific Limitations

### Higher Error Rate

Oxford Nanopore reads have a higher per-base error rate than PacBio HiFi, even with Q20+ chemistry. This can lead to:

- Lower variant calling confidence (QUAL scores)
- Increased false-positive frameshift calls near homopolymer runs
- Reduced sensitivity for small (1bp) insertions/deletions

### Minimap2 Preset

ONT data uses the `lr:hq` minimap2 preset (requires minimap2 >= 2.26). This preset is optimized for high-quality long reads but may produce different alignment characteristics than `map-hifi`, particularly at VNTR boundaries.

### Clair3 Model Selection

Clair3 automatically selects the appropriate model when `--platform=ont` is specified. ONT models may have different sensitivity/specificity trade-offs compared to HiFi models, especially for frameshift mutations in repetitive regions.

### Platform-specific validation

ONT development evidence includes three original simulated samples. That small
panel cannot establish general sensitivity or specificity. Report ONT separately
from HiFi, including normal controls, failures, ambiguity and uncalled alleles;
do not transfer HiFi accuracy estimates to ONT.

### NEGATIVE from R9-model ONT runs is not validated

The v0.16.1 caller-stage discordance gate (`calling.stage_discordance_min_af`,
`calling.stage_discordance_min_depth`) turns most Clair3 pileup-vs-applied
disagreements into INCONCLUSIVE instead of NEGATIVE, but it only covers
frameshifts that at least the pileup stage called. On older simulated ONT
amplicon reads called with the R9 Clair3 model (benchmark set `ms_ont_sub`),
one simulated pathogenic sample (`pair_5178`) still reports
NO_PATHOGENIC_VARIANT_DETECTED: its delGCCCA event is missed by every Clair3
stage, including the pileup stage, so there is no pileup-only call for the
gate to compare against. This residual false negative is a known, accepted
limitation of this release, not a regression; a caller-side remedy for events
that no Clair3 stage reports is deferred to a later engine change.

Two further pre-existing false negatives are unchanged by v0.16.1 and remain
open: `ms_pacbio_sub pair_5190` and the same simulated design reproduced on
both `leg_ont` and `leg200_ont` (`dev-clean-ont_amplicon_r10-0016`). In all
three cases both alleles are stage-concordant, so the v0.16.1 gate has no
discordant pileup evidence to act on; these are not caller-stage discordance
failures and are listed here as open, unresolved false negatives.

## Evidence status and reference fill

Equal lengths are not evidence of identical sequences. Unphased multiple variants,
missing genotypes and disconnected phase sets remain unresolved. An empty VCF
means no retained variants, not a confidently normal complete genotype. A single
observed length may still represent one or two biological alleles. Candidate
aliases in legacy allele fields are labeled and are not independent observations.

Consensus retains reference sequence at uncalled sites. Whole-VNTR read coverage
and uncovered-base masking are not yet validated; `reference_confidence` remains
`unverified`. Exact agreement with simulated truth does not prove that every base
had independent read support. Low-depth and partial reads require particular
caution when interpreting candidate sequences.

In v0.12.0, reference ladder generation and consensus trimming use the proximal
left flank (`flanking_left[-flank_length:]`) matching the primer-adjacent sequence
present in actual amplicon reads. Trimming derives its anchor from this proximal
sequence and records `exact_anchor` or fixed fallback status, ensuring exact anchor
alignment across contigs.

Length candidates still use secondary alignments for reference fit. `reads` is
retained as a legacy alias for `alignment_records`, not independent molecules.
`primary_alignment_records` counts records without QNAME deduplication, since
simulated distinct reads can share names. `molecule_count` remains unknown without
source/UMI evidence. Error-tolerant anchor diagnostics do not replace inference:
adjacent lengths, minority alleles and ONT anchor failures remain limitations.

## Experimental read-backed phase

For a single observed length with unresolved multiple heterozygous variants,
the Python `read_phase=True` option or JSON `calling.read_phase: true` can use
installed WhatsHap to link retained variants using remapped
primary records. Each record receives a temporary unique name; the source map
preserves original identities and does not deduplicate colliding names. Linked
biallelic SNPs and indels are tested with WhatsHap 1.7. Missing optional WhatsHap
is reported as unavailable; a present tool's failure propagates. No whole-sequence
reference confidence or molecule independence is inferred from a phase set.

Disconnected evidence, incomplete genotypes and unsupported linked multiallelic
sites remain unresolved. A false upstream length split bypasses same-length
phasing. Identical trimmed VNTR sequences do not inherit independent haplotype
credit from heterozygous flanking sites. This method cannot establish homozygosity
from absence of variants or repair unobserved reference-filled sequence.

This option is disabled by default and has no dedicated CLI flag: cached development
validation gained individual exact sequences but introduced an extra supported
mutation on a collapsed length cluster. A targeted remedy is under investigation.

The phasing read list is internally selected and may contain zero records for
one genotype haplotype. That haplotype is a genotype complement, not a sequence
observed on a selected read. A phase set does not establish two observed biological
alleles. `phasing_selected_records_by_phase_set` exposes these counts separately
from all primary input records; no validated minimum-support policy is applied.

`length_selection_evidence` lists subthreshold contigs and passing clusters not
selected by the current inference. Its counts measure alignment records, not
molecules. Null primary counts mean the BAM was not queried (including a missing
path); zero means no excluded mapped primary records. The summed subthreshold
count mixes distribution tails with possible separate length groups and must
not be interpreted as evidence for an additional allele without the per-contig
distribution and original-read assignment.

## Variant support and confidence

`vcf_support_status=projection_unavailable` means selected VCF alleles could not
be projected with a verified replay, including mixed IUPAC heterozygous indels.
It does not mean the variant was absent. The legacy boolean remains false until
exact concordance is proven; the HTML report labels this support unavailable.
Unavailable or ambiguous projection adds no support penalty to dictionary-fit
confidence, whereas actually absent support receives that penalty.
`vcf_support_status=heterozygous_genotype_unresolved` marks an event whose VCF
record is heterozygous and was replayed under `-H I`; it is never support and
receives no absence penalty. The separate
existing heuristic penalty near allele boundaries still applies.
These scores remain heuristic weights, not probabilities or empirical base support.


Custom reference layouts control constructed fixed repeats and trimming anchors.
Backward classifier recovery still uses dictionary category order before its
regular matching loop; arbitrary noncanonical layout reconstruction is not
scientifically validated. An explicit FASTA must match its generating dictionary
and flank configuration; the pipeline does not yet prove complete reference
compatibility from the FASTA alone.

## Hybrid Engine (Experimental)

`--engine hybrid` is **experimental and not the default engine**; `ladder`
remains the default until the benchmark decision rule is met on the sealed
test split. Every `hybrid.*` default (see the
[configuration guide](../guides/configuration.md#hybrid-engine-experimental))
is provisional and tuned on development/validation splits only. The numbers
below come from internal regression runs on a development branch (frozen
simulated panels, PRJEB92208, and in-house genomic data whose outputs stay
outside the repository); they are not part of the automated test suite and
are not a calibration or release claim.

### Detection limits

- **Equal-length heterozygotes near the `het_af_min` floor.** A minor allele
  fraction below `hybrid.het_af_min` (default 0.2, valid range [0.01, 0.5])
  never forms a phase-split candidate site. At the default, a minor allele at
  roughly 15-20% of an equal-length pair produces a silent `none` split (no
  flag), because `het_af_min` is set above `het_min_group` (0.15) by design.
  Detecting or at least flagging that edge needs either a lower `het_af_min`
  (calibration) or an explicit low-AF flag rule.
- **Equal-length heterozygotes with a single differing event.** When both
  alleles have the same repeat count and differ at one event only, the
  length model finds one peak and the linked-site split needs
  `min_linked_sites` (2) events. With the default
  `hybrid.phase_single_event_split = "indel"` a length-changing event (every
  frameshift, e.g. dupA/dupC) splits the peak into two alleles and the event
  is scored as usual. The reads of each allele are selected by the event
  itself, so the event's read support is conditional on that split; the
  independent evidence is the candidate-site test (run background, strand
  bias, `het_min_group`) and the peak-level share gate: the one-sided lower
  confidence bound (`phase_single_event_alpha`) of the stutter-deconvolved
  minor share, over a fresh seeded sample of at most
  `phase_single_event_bound_reads` reads, must reach `het_af_min`. On synthetic wild-type ONT-like reads (1200 reads) with a
  site-specific one-base excess at one C7 run on both strands, excess shares
  of 0.10-0.30 were never PATHOGENIC. From about 0.35 such an artefact cannot
  be told apart from a real minor allele of that share and is called like one.
  That stress test used a C7 run inside a canonical `X` unit; a symmetric
  excess of 0.35 or more at a C7 run in a variant unit or at the edge units of
  the array is untested. A
  single substitution-only event (not split at the default) and any split
  that yields identical alleles keep the phase basis `unconfirmed_single_site`
  (selection status `unresolved_single_site`) with the
  located reason `unresolved heterozygous site at repeat N`, which makes the
  result INCONCLUSIVE instead of NEGATIVE.
- **True dupC with high deletion stutter.** `hybrid.event_max_alternative_frac`
  (default 0.25) rejects a homopolymer event whose no-event mixture weight
  exceeds it. The mixture convolves each allele with the stutter profile of its
  own run length (`hybrid.hp_stutter_model = "length"`), so the heavier deletion
  stutter of a longer run is not taken for a no-event allele. MUC1 has no C8 run
  besides a dupC run, so the C8 profile is always extrapolated from the sample's
  shorter C runs (per-base growth capped at `hp_stutter_max_growth`). Where the
  sample's stutter saturates (real ONT "+" reads) the extrapolated C8 allele
  would be read as C7 about as often as C8; such a strand falls back to the
  shifted background (`hp_stutter_max_event_confusion`). On
  simulated HiFi reads with about 31-32% C7 reads at the true C8 run against
  7-8% deletion stutter at the C7 runs (a steeper step than the trend below C7),
  the estimated alternative share of the true dupC fell from 0.263/0.284 (the
  former shifted background) to 0.140/0.187: supported at the default, but
  closer to the limit than on real ONT amplicon reads (PRJEB92208 dupC
  positives, 0.07-0.12). The one consensus-error dupC in the simulated panels
  (a spurious C8 over a 52:47 C8/C7 read split) stayed `discordant` at 0.379
  (0.439 before), so the margin between true events and that error is 0.19
  (0.155 before). Tolerated wild-type share at the default limit (synthetic,
  three seeds, share of seeds still `supported`): log-linear HiFi-like stutter
  1/3 at 20% and none from 25%; the D1-shaped step none from 20%; the
  saturating ONT "+" shape (shift background after the guard) 3/3 up to 25%
  and none at 30%. A wild-type share of 25% or less on such ONT data can
  therefore pass as a pure dupC, with either background. A stutter step at
  the event length that is steeper than any trend in the sample's shorter
  runs cannot be predicted from those runs, and run length alone cannot tell it apart from a C7/C8 mixture;
  `event_max_alternative_frac` stays a calibration trade-off.
- **Residual single-base HiFi consensus misses.** On a 40-case frozen
  simulated panel (`simpanel`), 77/80 alleles were sequence-exact at commit
  `ca81a97` (see "Validation numbers" below); every mutation event and every
  allele count was still correct, and each inexact allele carried its own
  `residual_sites` flag. The misses are single-base indels in
  homopolymer-adjacent contexts
  (`GCG CCC G CA`, `GCG C5 A`) where the simulated reads' own majority differs
  from the simulator's ground truth by one base (a data property, not a
  decision error), plus one case where a substituted base and its neighbouring
  insertion slot are split across two separate pileup-vote columns. Neither
  failure mode changed a clinical call.
- **PRJEB92208 ONT amplicon negatives.** Most amplicon runs on that dataset
  produce several length peaks below `hybrid.min_peak_reads`/
  `far_peak_min_frac`/`near_peak_min_frac` from PCR-smear reads, so sample
  `selection_status` is usually `unresolved_rejected_peak`. That correctly
  blocks a NEGATIVE result (no false reassurance) but also means most
  amplicon runs without a pathogenic event end up INCONCLUSIVE rather than
  NEGATIVE: internally, about 6 of 9 amplicon runs reached a callable
  (non-`unresolved_rejected_peak`) primary-amplicon selection (measured at
  commit `2b0072b`, before the Task 13b fix; that fix (`git diff
  2b0072b..ca81a97`) touched only `hybrid/evidence.py`,
  `hybrid/polish.py` (the homopolymer consensus vote) and
  `settings_hybrid.py` (two new settings) -- not `hybrid/lengths.py`
  (length-peak selection) or `hybrid/allele_fields.py` (`selection_status`)
  -- but this figure itself was not re-run at `ca81a97`). No pathogenic call
  was produced on a known-negative sample.

### Validation numbers

Two internal validation passes describe this engine, both on branch
`feat/hybrid-engine`: one at commit `2b0072b` (before the Task 13b fix) and
one after it, at commit `ca81a97`. The Task 13b fix (`git diff
2b0072b..ca81a97`) changed `hybrid/evidence.py` (per-event read-level
support: the event vs. no-event vs. read-derived-alternative comparison and
the homopolymer stutter mixture fit), `hybrid/polish.py` (the homopolymer
consensus vote counts only reads that observe a run cleanly, so a read that
merges the run with its neighbour no longer votes) and `settings_hybrid.py`
(`event_context_units`, `event_max_alternative_frac`); it did not touch
`hybrid/lengths.py` (length-peak selection) or `hybrid/allele_fields.py`
(sample `selection_status`). The frozen simulated panels and the 9-library
PRJEB92208 **amplicon** run were re-executed afterward and are reported at
`ca81a97` below. The PRJEB92208 **genomic/WGS** invocations, the in-house
genomic runs, and the `clinical_benchmark.py`-scored sequence-exactness
check were run once, at `2b0072b`, and were **not** re-run after the fix --
each row below states which commit it was measured on.

**Re-verified after the fix (commit `ca81a97`):**

| Check | Result |
| --- | --- |
| Frozen simulated dev panel, alleles sequence-exact | 77/80 (mutations detected 24/24). One case names an extra, read-support-`discordant` variant alongside the correct causative event from a residual 1-base consensus error; the clinical decision (PATHOGENIC on the true event) is unaffected, but the panel's by-name mutation scorer counts it as a false positive. |
| Frozen simulated held-out panel, alleles sequence-exact | 77/80 (mutations detected 28/28) |
| PRJEB92208 amplicon dupC-positive controls (4 libraries) | PATHOGENIC with `read_support.status == "supported"` on all of them |
| PRJEB92208 amplicon non-dupC samples (5 libraries) | none reached PATHOGENIC (all INCONCLUSIVE) |

**Measured once, before the fix (commit `2b0072b`, Task 13) -- not re-run
afterward:**

| Check | Result |
| --- | --- |
| PRJEB92208 HG002 vs. an independent full-sequence assembly (`clinical_benchmark.py`-scored, the amplicon library and the genomic WGS invocation) | both literal sequence-exact |
| PRJEB92208 genomic/WGS invocations (2 libraries) | both INCONCLUSIVE, on insufficient/low spanning depth |
| In-house ONT genomic samples (5 libraries, local only) | no phantom fragment alleles; a low-spanning-depth dupC sample reported INCONCLUSIVE (insufficient depth), not PATHOGENIC or NEGATIVE |

These are development-branch regression runs, not a release validation: the
frozen-panel gate for this engine (>=80/80 sequence-exact) was not met, the
public PRJEB92208 benchmark record was not refreshed after the Task 13b
evidence fix, several checks above were never re-run after that fix, and the
sealed test split was never used to tune a default.
