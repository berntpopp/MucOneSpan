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

## Variant Calling

### Long VNTR Alleles (>100 repeats)

Clair3 has reduced sensitivity for detecting 1bp insertions within very long tandem repeat regions (>6kb). For alleles with >100 repeat units, single-nucleotide frameshifts (e.g., dupC) may not be detected by Clair3, leading to false negatives.

**Affected test samples:**
- `sample_dupc_100_120`: dupC mutation in 100-repeat allele not detected
- `sample_long_120_140`: mutations in 120+ repeat alleles not detected

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
