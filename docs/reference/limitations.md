# Known Limitations

!!! note "Engines"
    Since 0.17.0 the [hybrid engine](#hybrid-engine) is the default. The
    sections from "Allele Length Detection" to "Variant support and
    confidence" describe the **deprecated ladder engine** (`--engine ladder`,
    minimap2 + Clair3 + bcftools), which stays available until a later release
    removes it (see the [migration guide](../guides/migration.md)). Its
    measured behaviour is unchanged by 0.17.0 except for two gates: a ladder
    run where one allele's depth is `not_assessed` while another's is assessed
    is now INCONCLUSIVE (0.16.1: NEGATIVE) and a frameshift on the
    `not_assessed` allele is no longer PATHOGENIC; any unknown depth status
    fails closed when a `depth_basis` is present; and the unresolved-selection
    reason no longer prints `secondary mode fraction None`.

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

## Hybrid Engine

Since 0.17.0 the hybrid engine is the **default** engine for amplicon and
genomic input. Every `hybrid.*` default (see the
[configuration guide](../guides/configuration.md#hybrid-engine)) was tuned on
the benchmark development split and confirmed on the validation split; the
sealed test split never informed a default. The numbers below come from
internal regression runs (the MucSim-Bench v4 splits, frozen simulated panels,
PRJEB92208, and in-house genomic data whose outputs stay outside the
repository); each states the commit it was measured at. They are not part of
the automated test suite.

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
- **Heterozygous homopolymer runs below the split floor.** A run-length site
  becomes a candidate only above `max(het_af_min, phase_run_bg_multiplier x
  background)`. Under heavy stutter a real heterozygous run can stay below it:
  a simulated HiFi equal-length heterozygous dupC (1000x) showed 42% C8 reads
  against 13.7% C8 at the peer C7 runs, under the x4 floor of 0.55, and was
  merged into one wild-type consensus and reported NEGATIVE. This was found on
  the benchmark validation split, not the development split. A single
  unsplit length peak is now tested against a lower safety floor
  (`hybrid.phase_run_safety_multiplier`, default 2.0); a run above it makes the
  result INCONCLUSIVE (`unresolved_run_site`, located), never PATHOGENIC, so such
  a carrier is not called. Costs and limits:
  - Wild-type runs with site-specific stutter at or above `het_af_min` are
    flagged the same way. In the v4 development panels the tier flags 3 of the
    4 equal-length (single-peak) wild-type HiFi samples (runs at 21-27% of reads
    at one length) and none of the 3 ONT ones.
  - The tier keeps `het_af_min` as its share floor, so a heterozygous run whose
    minor length is seen in fewer than `het_af_min` of the reads (strong
    allele imbalance combined with heavy stutter) is still not flagged.
  - With two length peaks the tier does not apply: each peak is one allele.
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
- **PCR dimers and smear between the alleles.** A head-to-tail PCR dimer is
  read as `allele_a + junction + allele_b`, a length peak near `L_a + L_b`
  (`2 x L` when both copies come from one allele). It is recognised only on
  structural evidence: each read must carry an internal motif-9 -> motif-1
  amplicon junction, and both of its parts must fall in an accepted allele's
  window. Length alone never makes a peak a dimer, so a real allele at twice
  another allele's length stays a gate-relevant rejected peak. Limits:
  - A dimer whose junction motifs carry more than `anchor_max_edits` edits,
    a trimer, or a dimer of an allele that was not accepted is not
    recognised and keeps the sample INCONCLUSIVE.
  - More dimer reads than `dimer_max_parent_frac` x the smaller parent's
    support are not treated as a minority artefact.
  - Head-to-head (inverted) dimers produce no extra length peak and are not
    counted.
  - The smear test between the alleles runs only when the shorter allele is
    the top peak and only re-judges `support_below_threshold` candidates. A
    real minor allele between the alleles (contamination, mosaicism) that is
    not a significant excess over the local smear is called smear, the same
    detection floor as the below-top smear test.
  - Measured at commit `e324fa3` on the v4 benchmark (INCONCLUSIVE share of
    all cases, before -> after): development `standard` 0.233 -> 0.122 and
    `clean` 0.189 -> 0.078; validation `standard` 0.233 -> 0.178, `clean`
    0.089 unchanged. Every changed case was a normal going from INCONCLUSIVE
    to NEGATIVE with its allele lengths correct; no false positive and no
    NEGATIVE on a pathogenic case. A simulated long allele with PCR dropout
    (102 units, 1% of reads) stays a gate-relevant peak.
- **Equal-length normals with strong single-run stutter.** The
  `unresolved_run_site` safety tier above cannot tell a wild-type run with
  strong site-specific stutter from a heterozygous run, so such normals are
  INCONCLUSIVE, not NEGATIVE (v4 normals of all three sets moved NEGATIVE ->
  INCONCLUSIVE when the tier was added: 2 of 77 on the development split and 1
  of 77 on the validation split). The ratio is fixed and has
  no depth term, so a low-depth ONT normal can be flagged too (one validation
  ONT genomic normal with 57 phase reads). Where the background falls back to
  shorter same-base runs it is under-estimated, which flags more, never less.
- **Whole-unit PCR slippage.** Clusters one or two repeat units below an
  accepted allele (PCR slippage, a few percent of that allele's reads) are a
  significant excess over the smear and stay gate-relevant rejected peaks,
  so the sample is INCONCLUSIVE. This is the remaining cause on PRJEB92208
  HG001-HG004 and on three v4 validation cases (clusters one unit below an
  allele).
- **Inter-allele smear near long alleles.** The background of the smear test
  between the alleles can be inflated next to long alleles (above about 100
  units), which lowers the chance that a real minor allele there is flagged.
- **PRJEB92208 ONT amplicon runs.** Most amplicon runs carry clusters of
  PCR-product reads below or above the alleles, so `selection_status` is often
  `unresolved_rejected_peak`. That blocks a NEGATIVE result (no false
  reassurance) but leaves the non-dupC runs INCONCLUSIVE. The PRJEB92208 runs
  hold no dimer peak (at most one dimer read per run), and their smear between
  the alleles is recognised. HG001-HG004 stay INCONCLUSIVE: their remaining
  gate-relevant peaks lie below the top allele and are a significant excess over
  the smear (HG001 and HG004: a cluster two units below the top allele with
  2.1% and 2.7% of its reads; HG002 and HG003: clusters far below it, some
  `smear_ambiguous`), or are clusters of 3-11 reads a few units above the longer
  allele. No pathogenic call was produced on a known-negative sample.

### Validation numbers

**MucSim-Bench v4** (simulated; development split used for tuning, validation
split for confirmation only; measured on the `e324fa3` snapshot, whose hybrid
defaults are the 0.17.0 defaults). PATHOGENIC is the share of pathogenic
cases called PATHOGENIC; INCONCLUSIVE is the share of all cases.

| Split | Set | PATHOGENIC | INCONCLUSIVE | False positives | NEGATIVE on a pathogenic case |
| --- | --- | --- | --- | --- | --- |
| dev | standard | 53/57 = 0.930 | 11/90 = 0.122 | 0/33 | 0 |
| dev | clean | 54/57 = 0.947 | 7/90 = 0.078 | 0/33 | 0 |
| dev | clean2 | 17/19 = 0.895 | 2/30 = 0.067 | 0/11 | 0 |
| val | standard | 52/57 = 0.912 | 16/90 = 0.178 | 0/33 | 0 |
| val | clean | 54/57 = 0.947 | 8/90 = 0.089 | 0/33 | 0 |
| val | clean2 | 17/19 = 0.895 | 2/30 = 0.067 | 0/11 | 0 |

`clean2` is one case short of the 0.90 PATHOGENIC target on both splits. The
sealed test split has not been run.

**Frozen simulated panels** (commit `e324fa3`, same hybrid defaults;
PATHOGENIC/INCONCLUSIVE/NEGATIVE counts):

| Panel | Normals | Pathogenic |
| --- | --- | --- |
| `simpanel` (40 cases) | 0/2/14 | 24/0/0 |
| `heldout` | 0/4/8 | 26/2/0 |
| `ms_ont_sub` (ONT) | 0/1/38 | 38/1/0 |

No false positive and no NEGATIVE on a pathogenic case. Sequence exactness on
`simpanel` was last measured at commit `ca81a97`: 77/80 alleles sequence-exact
(the residual single-base misses described above); it was not re-measured at
the 0.17.0 defaults.

**PRJEB92208** (public ONT data; `benchmarks/clinical/prjeb92208/hybrid-engine.json`,
re-run with the 0.17.0 release candidate at commit `a185ecc`, 2 threads):

| Check | Result |
| --- | --- |
| Amplicon dupC-positive controls MP1-MP4 | PATHOGENIC, dupC `read_support` `supported` on all four (242/336, 396/652, 3899/6551, 304/512) |
| Amplicon HG001-HG004 | none PATHOGENIC (all INCONCLUSIVE, `unresolved_rejected_peak`) |
| Amplicon MP5 | INCONCLUSIVE (`unresolved_max_alleles`) |
| HG002 amplicon vs. an independent full-sequence assembly | both alleles literal sequence-exact |
| Callability (amplicon invocation, 11 runs incl. the 2 WGS runs) | 7/11 |
| Genomic/WGS invocations (`--assay genomic`, 2 runs) | both INCONCLUSIVE on spanning depth (HG002 28/15 spanning reads, `low`; the identity-unresolved MP1 WGS run 8/8, `insufficient`); HG002 alleles literal sequence-exact |

**In-house ONT genomic samples** (5 libraries, local only; measured once at
commit `2b0072b`, **not re-run** at the 0.17.0 defaults): no phantom fragment
alleles; a low-spanning-depth dupC sample was reported INCONCLUSIVE
(insufficient depth), not PATHOGENIC or NEGATIVE.

These are internal regression runs, not a release validation on an
independent cohort: reported-control labels do not establish independent
positive truth, healthy-sample labels do not establish endpoint-specific
negative truth, and the genomic numbers above cover few libraries.
