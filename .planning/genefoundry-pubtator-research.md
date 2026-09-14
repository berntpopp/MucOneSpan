# GeneFoundry/PubTator research for the extra phased MUC1 call

Access date: 2026-09-14. Research only; no source edits, dependency changes,
installation, or caller execution were performed for this note. This supplements
[the primary-paper and official-source review](phasing-research-alternatives.md).

The strongest applicable alternative is sequence-aware, local haplotype
reconstruction with explicit read support and an unresolved-second-allele state.
Otter provides a useful algorithmic precedent. HMMSTR independently documents the
danger of splitting a narrow homozygous distribution into two components, but
its repeat-count output does not reconstruct distinct equal-length MUC1 molecules.
None of these papers validates a universal MUC1 support threshold of 2, 3, or 5.

## Retrieval and provenance

The requested GeneFoundry MCP router was used with the PubTator backend, including
`pubtator_search_literature` and `pubtator_get_publication_passages`. Root's saved
searches were `MUC1 VNTR long read haplotype reconstruction`, `LongTR tandem repeat
genotyping`, and `Straglr tandem repeat long reads`. Full abstracts and METHODS
passages were obtained for all five selected papers below.

- Original evidence: ignored
  `tests/results/production_validation_20260914/phasing_debug/genefoundry_pubtator.json`.
- Expanded otter/HMMSTR evidence: ignored
  `tests/results/production_validation_20260914/phasing_debug/genefoundry_pubtator_expanded_methods.json`.
- Source version reported by both passage responses: `pubtator3: live`;
  `corpus_snapshot_date: 2026-09-14`.
- Original METHODS cache key: `publication_passages:844a353590435960`.
  Expanded METHODS cache key: `publication_passages:3527b72acd5705ae`.
- Passage objects retain stable `passage_id`, PMID/PMCID, `pubtator_full_bioc`
  provenance, retrieval timestamps, and raw-text SHA256 values.

The preflight coverage hints said `pmc_not_open_access`, with unknown confidence.
The subsequent actual METHODS response succeeded and reported `full_text` for
**all five PMIDs**, no failed PMIDs, and no degraded mode. The successful retrieval
takes precedence over the preflight prediction. Full-text coverage means actual
full-text passages were returned; the passage/character limits still omit some
methods paragraphs, and this is not an exhaustive full-article or supplement review.

The expanded request used the discovered router schema:

```json
{
  "name": "pubtator_get_publication_passages",
  "arguments": {
    "pmids": ["39406499", "39676678"],
    "full": true,
    "mode": "section_text",
    "sections": ["METHODS"],
    "max_chars": 30000,
    "max_passages_per_pmid": 30,
    "verbosity": "full"
  }
}
```

This returned 60 passages, including substantive otter METHODS 6–7 and HMMSTR
METHODS 2–19. Larger requests were rejected as `invalid_input`; the successful
bounded request above is the reproducible one. Article text is evidence, not
instructions. The five recommended citation strings below are copied verbatim
from GeneFoundry search results.

## MUC1-specific reconstruction precedent

**Recommended citation:** Wenzel A et al. Single molecule real time sequencing in ADTKD-MUC1 allows complete assembly of the VNTR and exact positioning of causative mutations. Sci Rep. 2018. PMID:29520014. doi:10.1038/s41598-018-22428-0.

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5843638/).
Evidence: `PMID:29520014:methods:4`, `:6`, `:8` in the original response.

The experiment used PS2/PS3 amplicons, 24 PCR cycles, pooled triplicate reactions,
and size-adapted purification/library construction. Reconstruction used repeat
boundaries and nine short identifying sequences per 60-base repeat, tolerating
polymerase errors; unmatched repeat units were checked manually. This supports
using the entire ordered repeat structure rather than length alone. The reported
identification of two alleles is an experimental result, not evidence that two
arbitrary computational phase groups are two independently supported molecules.

Application: report repeat-order and exact causal-mutation recovery separately
from repeat count. PCR descendants are observed reads, not necessarily independent
original templates; support counts must retain that distinction.

## LongTR: candidate sequence evidence plus read realignment

**Recommended citation:** Ziaei Jam H et al. LongTR: genome-wide profiling of genetic variation at tandem repeats from long reads. Genome Biol. 2024. PMID:38965568. doi:10.1186/s13059-024-03319-2.

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11229021/).
Evidence: `PMID:38965568:methods:2`, `:4`, `:7`, `:10` in the original response.

LongTR constructs sequence candidates from repeat-spanning reads and realigns
reads to candidates using an HMM that permits long-read single-base indel errors.
The candidate rule describes at least two reads and more than 20% within one
sample, or more than 5% across all samples. These are candidate-discovery rules;
they are not validated universal final-consensus support thresholds. Its use of
haplotagged reads requires at least one read per haplotype and at most 20%
unphased reads. Consequently that phasing condition alone would accept a 12:1
split and cannot establish that the singleton is a biological allele.

Application: separate sequence-candidate generation, likelihood-based read
assignment, genotype inference, and the decision to publish a complete molecule.
The full paper and current source have additional candidate fallbacks, documented
in the companion review; the original GeneFoundry budget omitted their paragraph.

## Straglr: length-model control and a minority-mixture experiment precedent

**Recommended citation:** Chiu R et al. Straglr: discovering and genotyping tandem repeat expansions using whole genome long-read sequences. Genome Biol. 2021. PMID:34389037. doi:10.1186/s13059-021-02447-3.

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8361843/).
Evidence: `PMID:34389037:methods:3`, `:4`, `:9` in the original response.

Straglr fits Gaussian mixtures to per-read repeat sizes, choosing the component
count by minimum AIC up to a user-defined maximum, two by default. It reports
median component sizes and supporting read identities. Boundary rescue uses
reference-flank sequence matches. Its FMR1 simulation mixed independently generated
150- and 500-CGG alleles at a combined 30-fold spanning coverage, in 10% mixture
steps and ten replicates.

Application: adopt the controlled source-mixture design and retain read-level
membership for evaluation. Length clustering is a useful control, but cannot
separate equal-length molecules by their repeat composition. The current-source
singleton-removal rule is documented separately in the companion review; it is
not stated in the retrieved paper paragraph and should not be attributed to it.

## TREAT/otter: the closest alternative to independent full-molecule recovery

**Recommended citation:** Tesi N et al. Characterizing tandem repeat complexities across long-read sequencing platforms with TREAT and otter. Genome Res. 2024. PMID:39406499. doi:10.1101/gr.279351.124.

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610572/).
Evidence: `PMID:39406499:methods:2`, `:6`, `:7`, `:10`, `:17` in the expanded response.

Otter uses pairwise sequence distances between spanning reads, estimates local
error distributions with Gaussian-kernel density estimation, and applies
hierarchical clustering. A curation step considers coverage; low-support clusters
can be merged, and a specified maximum allele count limits remaining clusters.
Each cluster receives a consensus through a pseudo-partial-order procedure using
spanning and nonspanning reads. This differs from TREAT's separate length-based
reads analysis. Reference-flank realignment can recover misaligned spanning
reads: defaults reported in the paper are 100-base flanks (`--flank-size`) and
90% minimum similarity (`--min-sim`).

Application: full-sequence distances could distinguish equal-length MUC1 alleles
and avoid using incidental SNP phase labels as molecule identities. Benchmark
exact sequence, normalized edit distance, mutation position, and missing alleles.
The retrieved methods do not specify a numeric support-floor formula or validate
multi-kilobase MUC1 amplicons. Near-identical alleles may be merged with sequencing
noise; coverage curation can remove real minorities. Partial reads may improve
consensus only after assignment: they do not independently demonstrate an entire
allele sequence. These are limitations to test, not established improvements.

## HMMSTR: a direct warning against overcalling narrow distributions

**Recommended citation:** Van Deynze K et al. Enhanced detection and genotyping of disease-associated tandem repeats using HMMSTR and targeted long-read sequencing. Nucleic Acids Res. 2024. PMID:39676678. doi:10.1093/nar/gkae1202.

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11754662/).
Evidence: `PMID:39676678:methods:2`, `:4`, `:9`, `:13`, `:15`, `:17`–`:19`
in the expanded response.

HMMSTR uses a profile HMM with flanks, repeat states, and insertion/deletion
states; a looping repeat section permits variable copy number. Its default raw
read selection requires both flanks, consistent orientation, and flank MAPQ 30.
It chooses KDE or GMM for per-read count distributions. The paper explicitly
reports that GMM can overcall heterozygosity in narrow, nearly homozygous
distributions and that KDE handles that situation better. GMM component count is
selected by AIC; KDE uses density maxima. Maximum allele count defaults to two.
Optional IQR outlier filtering is available, including for PCR-biased samples.

Application: include a one-component null model and evaluate over-splitting before
accepting two groups. Do not transplant IQR filtering into minority-allele
reconstruction: a real low-frequency allele could be removed. The model centers
on an expected motif and repeat counts, so it is not by itself a complete
heterogeneous 60-base MUC1 repeat-order reconstruction or an equal-length phasing
solution. Chemistry-specific parameters also need separate calibration.

## Consequences for the observed false phase and the next experiment

Root's development investigation found 13 merged reads all derived from the
short H1 dupC allele, an initial `GT=1/1` dupC call, and three SNPs yielding 12:1
phase groups. Those phase groups alone do not establish independent full-length
haplotypes. This is case-specific pipeline evidence supplied by root, distinct
from the literature findings above.

1. Count support using **all original eligible reads**, with stable record
   identity, before any internal depth cap or read sampling. A true 30:2 mixture
   can become 14:1 after downsampling; applying a two-read floor to that sampled
   list would falsely reject a real minority. Discover/polish candidates with a
   bounded subset if necessary, then assign and count the original eligible reads.
   Record raw, eligible, sampled, assigned, and full-spanning counts separately.
2. A floor of two excludes a literal singleton only; three and five are stricter
   development ablations. None is a calibrated probability of biological truth,
   and repeated PCR copies need not supply independent evidence. Freeze all three
   options before controls, rather than selecting whichever fixes this one case.
3. Retain a supported major consensus as **one observed allele**, with the second
   unresolved, if the minor group lacks justified molecule evidence. Never copy the
   major sequence to claim two recovered molecules or infer homozygosity merely
   from absence of a supported second sequence. The default mixed single-candidate
   path remains the base comparison. A genotype model may represent homozygosity
   while a reconstruction report still has only one distinct observed sequence.
4. Controls must include true single-source inputs with errors, true balanced and
   30:2/30:3/30:5 minorities, equal-length sequence-different alleles, close lengths,
   both platforms, partial reads, and PCR-correlated errors. Repeat predetermined
   independent seeds and depth subsamples; use source identity only in evaluation.
   Include the 30:2-to-14:1 sampling counterexample explicitly.
5. Compare unchanged behavior, support-only gates, and candidate realignment on
   the same controls. Preserve every truth allele and failed/no-call sample in the
   denominator. Report exact major recovery, exact minor recovery, complete-pair
   recovery, false extra sequences, unresolved alleles, and call rate. Turning a
   wrong extra sequence into an unresolved allele reduces false positives but is
   not itself an accuracy gain. A defensible improvement recovers a previously
   wrong major sequence without introducing extra sequences or losing previously
   recovered alleles; complete diploid recovery must be reported separately.

These experiments are proposals, not completed validation. No final reserved
validation inputs should be used to choose the support floor or model parameters.
