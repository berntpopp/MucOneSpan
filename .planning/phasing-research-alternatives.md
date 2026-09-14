# Tandem-repeat reconstruction alternatives for the phasing extra-call investigation

Accessed **2026-09-14**. Bounded primary-source review, not a systematic review.
Sources: original papers, official repositories/documentation and official GitHub
API version metadata. No installation, dependency change, source change, or caller
execution. Recommendations below are experiments, not validated MUC1 defaults.
WhatsHap/Clair3-specific review is handled separately by the coordinator.

## Main recommendation

Prioritize a **read-to-candidate sequence assignment audit followed by within-group
consensus**, retaining an explicit one-observed-sequence outcome. Use LongTR as
a sequence-reconstruction comparator; use Straglr as a length-only control.
Evaluate TRGT's actual pinned implementation before treating its `cluster` option
as sequence-based support for huge equal-length MUC1 alleles. A literature label
or two allele assignments does not prove independent haplotypes.

The strongest newly verified limitation is in released TRGT 5.1.0: its clustering
distance becomes length-only for long read pairs. Consequently it is not an
automatic remedy for this project's central equal-length/close-length failure.
This is verified source behavior; the effect on MucOneSpan is a hypothesis until
paired development runs measure it.

## Version and source record

| Tool | Version actually checked | Immutable source / original paper |
| --- | --- | --- |
| TRGT | v5.1.0, release 2026-06-10; main matches tag | [Release](https://github.com/PacificBiosciences/trgt/releases/tag/v5.1.0), commit `b21c6217eeaef9a080205932fe29f730dc00e31f`; [Dolzhenko et al., Nature Biotechnology, 2024](https://doi.org/10.1038/s41587-023-02057-3) |
| LongTR | README advertises v1.2; current source commit dated 2025-09-23; paper evaluated v1.0 | [Pinned repository](https://github.com/gymrek-lab/LongTR/tree/fd73e0d45f239601a12854275ab4efa588010caf); [Ziaei Jam et al., Genome Biology 25:176, 2024](https://doi.org/10.1186/s13059-024-03319-2) |
| Straglr | Source version 1.5.6; current commit dated 2026-03-25 | [Version file](https://github.com/BirolLab/straglr/blob/014f80bd44bc08c20ee63e91aac82fb61676d491/src/version.py); [Chiu et al., Genome Biology 22:224, 2021](https://doi.org/10.1186/s13059-021-02447-3) |

GitHub API endpoints checked: each repository's `/commits?per_page=1`; TRGT's
`/releases/latest` and `/tags?per_page=5`. The TRGT release tag and source commit
agree. Current source was checked because historical paper algorithms and current
implementations differ. No claim that a different binary already installed here
has these versions; none was installed or executed in this research task.

## TRGT: useful read outputs, but concrete MUC1 hazards

The published approach extracts repeat-spanning reads, clusters pairwise edit
distances using Ward linkage, and produces allele consensus sequences. That is
the useful architectural idea: cluster actual observed repeat sequences and keep
read-to-allele assignments alongside consensus. The paper describes HiFi inputs;
it is not evidence of ONT applicability. [Original TRGT methods](https://pmc.ncbi.nlm.nih.gov/articles/11921810/)

**Released source caveat:** `get_dist` uses absolute length difference rather than
edit distance whenever `len(seq1)*len(seq2) > 10000`. Multi-kilobase MUC1 read
pairs necessarily reach that branch; same-length differences are invisible to
this distance. Additionally, consensus-length separation below 100 bp combined
with support imbalance greater than 4:1 triggers reassignment by alternating
read index and rebuilding two consensuses. This can affect a minority allele one
60-bp repeat away. The no-split fallback also alternates reads; one diploid read
can yield two copies. Cluster labels therefore cannot certify distinct molecular
sources. These are direct observations from v5.1.0, not inferred from its README.
[Pinned clustering source](https://github.com/PacificBiosciences/trgt/blob/b21c6217eeaef9a080205932fe29f730dc00e31f/src/trgt/genotype/genotype_cluster.rs)

Documented targeted settings: `--preset targeted` selects `--genotyper cluster`,
`--flank-len 200`, `--max-depth 10000`, flank identity fraction 0.8, alignment
scoring `1,0,1`, and disables the normal `rq>=0.98` filter. Default WGS mode uses
size genotyping, 250-bp flanks and depth 250. Thus the targeted preset bundles
several interventions; it is unsuitable as an isolated clustering ablation.
[Pinned CLI defaults](https://github.com/PacificBiosciences/trgt/blob/b21c6217eeaef9a080205932fe29f730dc00e31f/docs/cli.md)

PureTarget documentation explains allele dropout when expanded molecules fail
HiFi quality requirements, and recommends including retained fail reads. It also
documents a non-disableable low-quality-read purity heuristic near 90% similarity.
For heterogeneous MUC1 repeat compositions, discarded-read counts need inspection.
This PCR-free/targeted workflow recommendation does not demonstrate validation
for the current PCR amplicon design; unavailable original fail/subreads cannot
be recreated from retained BAMs. [PureTarget workflow and caveats](https://github.com/PacificBiosciences/trgt/blob/b21c6217eeaef9a080205932fe29f730dc00e31f/docs/puretarget.md)

If used as a comparator, retain `.spanning.bam`: it contains used reads,
`AL` assignments, locus identity and flank offsets. `AL` labels are tool allele
indexes, not simulator H1/H2. The `deepdive` command realigns used reads to allele
consensuses for visual support inspection. These derived alignments are not
independent validation of the sequences that generated them.
[BAM tags](https://github.com/PacificBiosciences/trgt/blob/b21c6217eeaef9a080205932fe29f730dc00e31f/docs/bam_files.md),
[Deepdive](https://github.com/PacificBiosciences/trgt/blob/b21c6217eeaef9a080205932fe29f730dc00e31f/docs/deepdive.md)

## LongTR: stronger sequence comparator with minority caveats

LongTR generates candidate sequences, including partial-order-alignment (POA)
consensus clusters, then scores candidate diploid pairs by read realignment with
a hidden Markov model. The paper's initial exact-sequence candidates require at
least two reads plus >20% within one sample or >5% across samples. If >25% of a
sample's reads lack a representative, fallback clustering starts at edit-distance
threshold 10, increases through 20/50/80/100/150/200/300/400/500/600/700, and seeks
>80% read representation. Clusters must exceed `min(10,0.1*n_excluded)` support.
Unphased genotype scoring assumes equal haplotype mixture weights; external
haplotags are used only with both haplotypes represented and ≤20% unphased reads.
Inference: minority groups can be missed before genotype scoring, and aggressive
distance merging can combine biologically close alleles. These mechanisms need
explicit testing rather than copying their thresholds. [LongTR methods](https://doi.org/10.1186/s13059-024-03319-2)

Operational defaults would exclude much of this panel: `--max-tr-len 1000`,
`--min-reads 10`, `--min-mean-qual 30`, `--min-mapq 20`. The catalog uses a 1-based
start despite being called BED. `--phased-bam` requires genuine haplotags; it does
not create them. `--indel-flank-len` defaults to 5. Capture `DP`, `DFLANKINDEL`,
`Q`, `PQ`, `GLDIFF`, `PDP`, and `INEXACT_ALLELE`; POA-derived versus directly
observed sequences remain distinguishable. [Pinned LongTR README](https://github.com/gymrek-lab/LongTR/blob/fd73e0d45f239601a12854275ab4efa588010caf/README.md)

A primary ONT benchmark used LongTR v1.2 with alignment parameters
`-1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1.0,-1.0`. It reports complete
dropout of its lower-quality disease cohort under default depth/quality filters,
and separately tested `--min-reads 1 --min-mean-qual 1`. This supports measuring
filter exclusions, not adopting one-read diploid calls as evidence. Its accuracy
results are from other loci/cohorts and cannot be transferred to MUC1.
[Primary ONT benchmarking study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13001349/)

## Straglr: independent size control, not full-sequence phasing

Straglr genotypes repeat lengths with GMM clustering. Current source selects by
AIC (not BIC), considers up to target-cluster-count plus four components, removes
or merges small/nearby clusters, and limits the retained cluster count. Thus
allowing two clusters is not itself proof that both were observed independently.
Length-only grouping cannot separate sequence-distinct alleles of equal length.
[Pinned clustering implementation](https://github.com/BirolLab/straglr/blob/014f80bd44bc08c20ee63e91aac82fb61676d491/src/cluster.py)

**README/source discrepancy:** current executable defaults include
`--min_cluster_size 0.1` (fraction), `--max_num_clusters 2`,
`--min_cluster_d 10`, `--max_check_size 5000`, and `--flank_size 80`.
Clusters with median size ≥5 kb bypass the minimum-cluster-size check. The
README's minimum-cluster-size default of two is stale for this pinned source.
Discovery defaults also cap motif length at 50 and require expansion ≥100 bp;
these are unsuitable discovery filters for 60-bp MUC1 repeats and a 60-bp
allele separation. Use an explicit locus catalog and `--genotype_in_size` for a
size comparator; do not mistake a rounded copy count for ordered repeat truth.
[Executable arguments](https://github.com/BirolLab/straglr/blob/014f80bd44bc08c20ee63e91aac82fb61676d491/straglr.py),
[Clustering behavior](https://github.com/BirolLab/straglr/blob/014f80bd44bc08c20ee63e91aac82fb61676d491/docs/clustering.md)

## Genotype phase and independent complete molecules are different contracts

Proposed project interpretation:

1. A phased genotype connects variant alleles within a supported phase block;
   it does not establish that every base of a reconstructed VNTR was observed.
2. A flank heterozygote can label chromosome origin if reads physically connect
   it to the repeat. It does not establish two different repeat sequences when
   their distinguishing in-repeat bases are absent or unsupported.
3. A read spanning both repeat boundaries gives full-molecule sequence evidence;
   partial reads give local evidence only. A consensus stitched from disconnected
   partial reads can contain an unsupported cis arrangement.
4. PCR copies are not independently sampled molecules merely because they have
   separate records. Existing duplicate QNAMEs also do not identify duplicate
   molecules. Preserve record IDs; use simulator source truth only in evaluation.
5. Two consensuses obtained by arbitrary splitting of a collapsed group can
   reproduce the same input error twice or introduce an extra mutation. Length
   clusters, haplotags, allele labels, and independent source evidence must be
   recorded as separate quantities.

These are proposed evidence contracts/inferences for this project; none of the
reviewed tools establishes all five automatically for PCR amplicons.

## Bounded experiment proposal; no parameter selection on final outcomes

**A. Diagnose the extra call before changing a caller threshold.** For the
specific false event, compare all source-assigned development reads against the
event-bearing candidate and its reverted-parent candidate at the actual locus.
Record full-repeat spanning status, per-base event coverage, alignment ambiguity,
strand, allele assignment and likelihood difference. Determine whether the event
comes from erroneous clustering, consensus stitching, a read error, reference
projection, or a real but wrongly assigned minority sequence. Keep calls and
no-calls in the same outcome table. This is an experiment proposal, not a claim
that a particular mechanism has already caused the extra call.

**B. Freeze locus geometry from reference/primer evidence.** All comparators must
use one explicit locus/reference coordinate system rather than interpreting the
multi-contig ladder as independent genomic loci. Validate both true flank
sequences and the actual primer-limited available span. Do not require 200–250
unique flanking bases if the amplicon does not contain them, or reduce an anchor
threshold until final calls improve. Freeze the physically supported anchor
length before seeing outcomes; report reads lost at each flank/quality step.
Keep total repeat count, variable block, and external tool locus boundaries
separate in the evaluator.

**C. Three fixed algorithm arms on development data.** (i) Current frozen caller;
(ii) LongTR unphased sequence reconstruction with an explicit `--max-tr-len 10000`
for this ≤140-repeat design envelope, default minimum depth/quality first, and
the published ONT alignment settings as a separately labeled ONT arm; (iii)
TRGT HiFi and Straglr size comparison, treated as diagnostics with the above
source limitations. Before any run, freeze whether a permissive quality/depth
diagnostic is included; do not replace the primary arm after seeing no-calls.
Report unfiltered outputs and filter-exclusion counts. There is no current
published parameter set that demonstrates these MUC1 amplicon requirements.

**D. A new reconstruction hypothesis, if comparator limits reproduce.** Extract
anchored sequences; propose medoid/POA candidates; realign every full read against
all candidates; compare one-sequence and two-sequence models. Treat assignment
ties as ambiguous rather than alternating them into artificial haplotypes. Use
both sequence evidence and robust length; do not replace edit distance by its
length lower bound for long reads. Fit an explicit allele-fraction parameter
instead of assuming 50:50 PCR yield. Validate complexity penalties and minimum
independent support on development normal/minority controls. These are proposed
design choices, not validated likelihood calibration or ready production code.

**E. Fixed challenge matrix and acceptance.** Include 60/60 identical, 60/60
in-repeat-different, 60/61 and 60/62, 25/140, each minority direction, cis/trans,
flank-only heterozygosity, disconnected phase sets, partial-only evidence,
duplicate QNAMEs, and low depth on each platform. Use existing exposed data for
parameter development only. When captured source labels exist, report per-allele
assignment precision/recall and unassigned fractions, not just majority-dominated
cluster purity. Evaluate exact sequence, ordered repeats, counts, event identity
and phase separately. A lost mutant subgroup must not be hidden by fewer calls.

Freeze one chosen configuration before untouched validation. Require no new
normal mutation-positive sample, no lost previously exact count/sequence/event,
and no reduction of call rate for an accuracy-improvement claim. Count unresolved
or failed samples and missing alleles in every appropriate denominator. A stricter
no-call policy can be useful but is a safety-contract change, not accuracy gain.
Do not tune on the reserved final panel; if it is exposed during debugging,
reclassify it as development and reserve a new panel.

## Addendum: coordinator's 13-read, 12:1 false split

The coordinator supplied this development observation after the source review:
all 13 merged reads originated from the short H1 dupC allele; the initial dupC
genotype was 1/1; three SNPs produced a 12:1 phase split, with the singleton
carrying errors. The tested phase parameters did not change that split. This
research task did not independently reproduce those pipeline runs.

**Source-backed minimum-support precedents, not a calibrated MUC1 rule:**

- Straglr removes singleton clusters regardless of its configurable fractional
  support threshold. Its prior singleton merge can assign a singleton to a
  neighboring cluster; this should not be interpreted as confirming its variant.
  [Singleton removal in source](https://github.com/BirolLab/straglr/blob/014f80bd44bc08c20ee63e91aac82fb61676d491/src/cluster.py)
- TRGT's candidate split requires the smaller group to have at least
  `max(2,round(0.01*n))` reads: two when n=13. The source's subsequent arbitrary
  splitting/assignment means this is a split-support precedent, not a validated
  rule for certifying two haplotypes.
  [Split criterion](https://github.com/PacificBiosciences/trgt/blob/b21c6217eeaef9a080205932fe29f730dc00e31f/src/trgt/genotype/genotype_cluster.rs)
- LongTR's exact-sequence candidate rule starts at two reads, but its use of
  external haplotags needs only one per haplotype. Therefore an existing LongTR
  phase tag is not sufficient evidence against the observed 12:1 failure.
  [Candidate versus phase prerequisites](https://doi.org/10.1186/s13059-024-03319-2)

**Proposed outcome:** retain the read-supported major candidate as one observed
sequence, with the minor candidate withheld from independently resolved allele
outputs. Keep its rejected sequence/evidence in diagnostics. State that the
second genomic allele is unresolved or unobserved; do not copy the major
sequence into a second allele, label the sample homozygous, or treat a singleton
phase block as independent haplotype proof. The majority is chosen by a fixed
read-evidence rule, never by knowing which candidate matches truth.

This is consistent with the evaluator's explicit unresolved-alias contract and
can preserve useful reconstruction instead of discarding both candidates. It
must retain the baseline mixed single-candidate path as a comparator and avoid
replacing it unless support and exact outcomes justify the change.

| Predeclared floor per independently resolved candidate | Meaning for 12:1 | What it does not establish |
| --- | --- | --- |
| 2 primary records | Rejects the singleton-only second sequence; leaves major eligible | Two coincident systematic/PCR errors may still create a false sequence |
| 3 primary records | Same result here; also withholds true two-read minorities | No universal probability of correctness without control data |
| 5 primary records | Same result here; also withholds true two-to-four-read minorities | Higher purity may simply reflect lost minority calls |

The floor is an **absolute read-support requirement**, distinct from minimum
locus depth or a relative allele-fraction threshold. A fourfold or 10% fraction
rule would confound PCR imbalance with lack of an allele. Do not count duplicate
alignments of the same record as additional support, and do not collapse distinct
records merely because QNAMEs collide. In PCR data, primary-record count remains
an observation count rather than a demonstrated count of independent molecules.

**Calibration experiment:** freeze floors 2/3/5 as three development arms with
all other settings identical. Test normal single-source controls that contain
one/two/four erroneous reads, and true second-source controls with one/two/three/
five/ten usable reads across equal/one-repeat-separated/well-separated lengths.
Include correlated-error/PCR-copy controls where source evidence exists. Apply
the same candidate-local support, full-span and phase-connectivity checks in
every arm. Report source assignment precision and minor-source recall separately,
major sequence exactness, false additional haplotypes, event FP/FN, whole-sample
exactness and missing-allele/call rates. Select using development controls only;
without such controls describe the floor as uncalibrated and keep it experimental.

An improvement from zero to one exactly recovered sequence out of two can be a
real major-allele reconstruction gain, even while the sample stays incomplete.
It is not whole-diploid recovery. Claim accuracy improvement only if the paired
comparison also preserves existing exact sequences/counts/events, introduces no
normal FP, and does not reduce the fixed call-rate endpoint. Report the unresolved
second allele and all losses explicitly. Withholding the singleton alone is a
safety-contract change; it does not establish an accuracy gain by itself.

## Limits of this review

The papers study genome-wide/targeted repeat collections rather than this exact
MUC1 PCR simulator and source-identifiability contract. Version-specific source
behavior matters; tool defaults and README descriptions can disagree. This review
did not run any tool, benchmark performance, assess a clinical cohort, or establish
the cause of the current extra call. It identifies testable alternatives and
several concrete reasons that a blanket threshold reduction or generic two-way
clustering may fail.
