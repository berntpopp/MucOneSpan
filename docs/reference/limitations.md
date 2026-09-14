# Known Limitations

## Allele Length Detection

### Homozygous Same-Length Alleles

When both alleles have the same repeat count (e.g., 60/60), the indel-valley splitting algorithm may incorrectly separate reads into two clusters. The pipeline reports `same_length: true` and uses a disambiguation strategy, but accuracy is reduced compared to heterozygous samples with distinct allele lengths.

### Extreme PCR Bias (Asymmetric Alleles)

For highly asymmetric allele pairs (e.g., 25/140), the shorter allele amplifies much more efficiently during PCR. The longer allele may have very low read coverage (<10x), making allele detection unreliable. The pipeline requires a minimum coverage threshold (default: 10 reads) per allele.

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

The bundled ladder currently uses the first 500 bases of the stored 10 kb left
flank, which is distal to the VNTR. This biological reference limitation is retained
to avoid an unvalidated reference change. Trimming now derives its anchor from the
same flank prefix actually used by the ladder and records exact-anchor or fixed
fallback status. This repairs the reference/trim mismatch, but does not establish
biologically correct left-flank reconstruction. Replacing the distal reference
requires a paired reference ablation and regenerated bundled resources.

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
confidence, whereas actually absent support receives that penalty. The separate
existing heuristic penalty near allele boundaries still applies.
These scores remain heuristic weights, not probabilities or empirical base support.


Custom reference layouts control constructed fixed repeats and trimming anchors.
Backward classifier recovery still uses dictionary category order before its
regular matching loop; arbitrary noncanonical layout reconstruction is not
scientifically validated. An explicit FASTA must match its generating dictionary
and flank configuration; the pipeline does not yet prove complete reference
compatibility from the FASTA alone.
