# Cached production-path read phasing validation — 2026-09-14

## Decision

The read-phasing helper passes its bounded software tests and improves exact
individual haplotype recovery on these development inputs, but **fails the
scientific no-regression gate**: a positive sample gains one extra supported
mutation annotation. Keep read-backed phase experimental and disabled by default;
do not advertise full default read-backed reconstruction as delivered. The root
coordinator is implementing that policy; this experiment did not tune thresholds.

The user subsequently requested investigation of whether this regression can be
fixed. Root-cause work remains open; no speculative fix is represented as proven.

## Experiment and provenance

Ignored reproducible script and complete artifacts:
`tests/results/production_validation_20260914/read_phasing/production_cached/`.
Command:

```
uv run --locked --all-extras python tests/results/production_validation_20260914/read_phasing/production_cached/run.py
```

Five pre-existing `same_length` cached cases were selected before scoring. Every
arm includes every expected sample. `provenance.json` records hashes of the same
filtered VCF, remapped BAM and actual contig reference supplied to both controlled
arms. Current development allele metadata supplies the same pre-existing inferred
lengths and contigs in both arms. Mapping, Clair3 and VCF filtering are reused,
not rerun. No truth is supplied to phasing, candidate selection or classification.

The controlled arms execute production `disambiguate_same_length_alleles` with
only the upstream mapping/calling/filter boundaries replaced by cached paths;
one arm disables the new phaser and the other runs it. Both execute actual
`build_consensus_per_allele`, classification, exact VCF sequence-concordance
validation, and strict `load_truth` / `load_observation` / `evaluate_sample` /
`aggregate`. Every emitted candidate gets full/trimmed FASTA, consensus context,
alleles, repeat JSON/text, summary and explicit cached-reconstruction run status.

The experiment used the initially active hook. Its script must explicitly enable
the experimental library option if replayed after the default-disabled policy
change. Existing artifacts record the actual behavior and remain unchanged.

Additional comparison arms load the historical saved caller outputs and the
saved pre-read-phasing development outputs. Old saved output has one invalid
artifact (`sample_dupc_60_80__n20_seed1702`: classification/allele identifiers
disagree); strict evaluation keeps its truth denominators and grants no recovery
credit. It was not repaired or dropped to make a favorable historical comparison.

## Controlled phase-only results

Five samples, ten truth haplotypes and two truth mutation events in both arms.
All metrics below have identical lower/upper assignment bounds.

| Endpoint | Phase disabled | Read phased |
| --- | ---: | ---: |
| Exact individual sequences | 1 / 10 | 4 / 10 |
| Exact individual repeat structures | 1 / 10 | 4 / 10 |
| Exact reported repeat counts | 4 / 10 | 4 / 10 |
| Exact complete diploid samples | 0 / 5 | 0 / 5 |
| Predicted alleles | 5 | 9 |
| Missing truth alleles | 5 | 1 |
| Extra predicted alleles | 0 | 0 |
| Exact annotation TP / FN / FP | 2 / 0 / 0 | 2 / 0 / 1 |
| Supported exact annotation TP / FN / FP | 2 / 0 / 0 | 2 / 0 / 1 |
| Exact annotation precision | 2 / 2 | 2 / 3 |
| Samples with any output / no output | 5 / 0 | 5 / 0 |
| Completed / ambiguous reconstruction | 0 / 5 | 4 / 1 |
| New normal mutation-positive samples | — | 0 |

“Completed” is the evaluator execution/evidence category, not proof of a correct
full diploid reconstruction. The three normal candidates emitted no mutation
alarms. This small panel does not establish population specificity or empirical
reference-base confidence. Missing alleles and extra mutation annotations are
separate endpoints: zero extra allele count does not erase the extra event.

| Sample | Exact sequence before → after | Mutation TP / FP after | After count errors |
| --- | --- | --- | --- |
| normal60_80 n20 seed1703 | 0 → 1 | 0 / 0 | 0, -20 |
| dupc60_80 n20 seed1701 | 0 → 1 | 1 / 1 | 0, -20 |
| normal60_80 n20 seed1701 | 0 → 1 | 0 / 0 | 0, -20 |
| dupc60_80 n20 seed1702 | 1 → 1 | 1 / 0 | 0; one missing allele |
| normal60_62 previously exposed heldout | 0 → 0 | 0 / 0 | +1, -1 |

The previously exact haplotype is retained, and mutation TP/FN and sample call
rate do not regress. Nevertheless the extra supported event is a material
regression and prevents unconditional accuracy promotion.

## Historical comparison

| Arm | Exact seq / structure / count | TP / FN / FP | Predicted / missing / extra |
| --- | --- | --- | --- |
| Historical saved | 0 / 0 / 3 | 1 / 1 / 1 | 8 / 2 / 0 |
| Saved current development | 1 / 1 / 4 | 2 / 0 / 0 | 5 / 5 / 0 |
| Controlled current phase disabled | 1 / 1 / 4 | 2 / 0 / 0 | 5 / 5 / 0 |
| Controlled current read phased | 4 / 4 / 4 | 2 / 0 / 1 | 9 / 1 / 0 |

Historical differences combine multiple production changes and one invalid
artifact; only the two controlled arms isolate read phasing. Saved development
and controlled-disabled metrics agree on this panel.

## Extra event requiring investigation

Both emitted candidates in `sample_dupc_60_80__n20_seed1701` contain
`X:dupC` at total-repeat index 25. Both receive `frameshift=true`,
`template_match=true`, `vcf_support=true`,
`vcf_support_status=exact_sequence_concordance`, QUAL27.24. Exact sequence
assignment pairs one length60 candidate to the length60 mutant truth and the
other length60 candidate to the length80 normal truth, where the duplicate event
is false. Their truth distance matrix is `[[0, 3], [1215, 1212]]`.

The current phase producer preserves VCF genotype allele identity; the possible
causal roles of initial collapsed allele assignment, genotype calls and phase
partitioning are being investigated separately. No phase or allele is selected
using mutation truth. Whole-sequence VCF concordance is not independent molecular
validation and does not prevent a wrong biological haplotype assignment.

## Software evidence and limits

Focused tests: 14 deterministic unit tests and 17 actual WhatsHap/bcftools
integration tests pass. Ruff and configured mypy pass for the new module.
Integration tests include exact cis/trans/shared SNP/indel sequences, collision
retention and source-pure assignments, name sensitivity, partial/disconnected
reads, missing GT/PS, empty variants and multiallelic limitations. Global CI and
packaging checks are coordinated by root.

Cached phasing uses no new sequencing, no Clair3 rerun and no unseen samples.
No cached read-haplotype source map is available. Exact allele/event scoring is
valid, but assignment purity cannot be inferred retrospectively from QNAME or
sequence similarity. Canonical60_60 false split remains a separate upstream
limitation and was excluded from this five-case experiment.

## Root-cause trace and parameter interventions

Executed diagnostic script `production_cached/debug_extra.py`; all raw records,
VCF calls, exact source joins and six parameter-trial outputs are retained under
`production_cached/extra_event_debug/`. No production code was changed.

1. The mutation exists **before phasing** at `contig_51:1992 G>GC`, QUAL27.24,
   genotype `1/1`. All 13 selected primary records directly support `GC` at this
   anchored CIGAR position. The mutation itself is not heterozygous or a phase
   switch; a genotype-preserving phaser must apply it to both GT selections.
2. The only phase-informative calls are SNPs `3420 C>A` (QUAL7.47), `3460 C>G`
   (QUAL22.88), `3462 C>A` (QUAL19.70), initially `0/1`. WhatsHap yields `1|0`
   at all three, common PS3420. Its read list assigns 12 selected records to
   GT1 and one record to GT2. The singleton is `S/66/ccs`, merged ordinal12,
   original-input ordinal14 (all ordinals in this paragraph are one-based).
3. Every merged record uniquely matches an original record by exact
   orientation-normalized sequence **and quality**, without using QNAME as a
   join key. The original ordinals are
   `[15,2,3,4,5,6,7,8,9,10,13,14,17]`. All 13 selected sequences are short,
   3606–3644 read bases. Thus QNAME collisions do not explain this extra event.
   Root independently observes omitted long-read evidence in the original input.
   These matches prove source-record retention, not simulator haplotype origin.
4. Both phased trimmed sequences contain 3601 bases. GT1 exactly equals the
   known mutant 60-repeat haplotype; GT2 differs from that sequence by only
   three substitutions, at trimmed one-based positions2921/2961/2963. Neither
   reconstructs the normal80-repeat truth. GT2 is assigned to that missing
   truth by the evaluator's fixed one-to-one sequence assignment, giving the
   additional dupC false annotation; it is not an extra allele-count error.
5. Actual trials with default coverage15, coverage10, coverage20, MAPQ5,
   `--distrust-genotypes`, and
   `--distrust-genotypes --include-homozygous` all retain exactly the same three
   `1|0` calls, PS3420 and every original genotype. No parameter tested fixes the
   missing normal80 reconstruction. Genotype-changing options remain diagnostic
   only and are rejected by the production helper's genotype-preservation guard.

The installed `phase --help` and version1.7 source explicitly support the two
genotype options. The current online guide's old “Trusting the variant caller”
paragraph conflicts with the project's release history, which describes genotype
likelihood-aware changes and the `--include-homozygous` option. For this
experiment, actual local help, source and executed outputs take precedence over
that stale paragraph. Sources: [official guide](https://whatshap.readthedocs.io/en/latest/guide.html),
[official release history](https://whatshap.readthedocs.io/en/latest/changes.html).

The causal failure is therefore not a swapped phase label or a mutation invented
by WhatsHap. The pipeline treats one selected short-read cluster as though it
necessarily contains both biological homologs, then promotes residual heterozygous
calls to two reconstructed biological alleles. Correcting input length/record
assignment and genotype/model evidence is required to recover the missing normal
allele. The exact contributions of sequencing error versus repeat-alignment error
to the three minor SNP calls are not established by these data.

Unsound shortcuts are excluded: deduplicating the same mutation across output
haplotypes loses genuine shared mutations; forcing homozygous ALT to heterozygous
contradicts both this selected read evidence and valid homozygous mutation cases;
choosing the mutation-free GT using truth is invalid. A support floor can preserve
a well-supported major candidate, but is a selective evidence rule, not recovery
of the missing long allele. That hypothesis is evaluated separately below.

## Support-floor experiment and decisive downsampling counterexample

Executed `production_cached/support_gate.py`. Full per-sample and per-fixture
results are in `production_cached/support_gate/`. This is a diagnostic selection
experiment only; it does not alter production outputs or delete candidates.
It uses WhatsHap read-list groups0/1 mapped explicitly to GT selectors1/2.

On the five cached cases, thresholds2,3,5 give identical results: four exact
sequences and structures, four exact counts, TP2/FN0/FP0, six predictions,
four missing alleles, zero extra alleles, four ambiguous reconstructions and
one completed reconstruction. All five samples retain an output. The three
removed outputs are singleton read-list groups. The normal60_62 candidate has
selected groups9:6 and retains both; the no-heterozygosity dupC case is unchanged.
Thus the majority exact sequences can be retained while avoiding the extra
annotation, but the normal80 alleles remain missing and three outputs are
explicitly withheld. This is not full diploid recovery.

Forty-five known-sequence real-tool fixtures evaluate true cis, trans, shared
homozygous insertion plus two heterozygous sites, shared insertion plus one
heterozygous site, and single multiallelic1/2. Each runs at input record ratios
10:1,10:2,10:3,10:5,1:1,30:1,30:2,30:3,30:5. All45 original genotype-selected
pairs are exactly correct. Each threshold is applied to the resulting tool read
list, producing135 explicit candidate-selection observations. The single-site
and single-heterozygosity shared fixtures bypass phasing and remain unchanged;
this is an important scope limitation, not evidence of adequate read support.

| Actual independent input records | Selected read-list records, all three multi-site designs |
| --- | --- |
| 10:1 | 10:1 |
| 10:2 | 10:2 |
| 10:3 | 10:3 |
| 10:5 | 10:5 |
| 1:1 | 1:1 |
| 30:1 | 14:1 |
| 30:2 | 14:1 |
| 30:3 | 13:2 |
| 30:5 | 13:2 |

The 30:2 case is decisive: threshold2 on WhatsHap's selected read list discards
a true minority haplotype that has two independent input records and was
previously reconstructed exactly. Threshold3 similarly discards the true
five-read minority in30:5. The read list is an algorithmically downsampled subset,
not a count of all compatible BAM records. A selection policy cannot substitute
these two denominators. Even without downsampling, a threshold deliberately
withholds correctly recovered true minority haplotypes below its floor;1:1 yields
no output under every tested floor. That tradeoff must not be hidden as an
accuracy improvement.

For true shared homozygous mutations plus two heterozygous sites, both adequately
supported sequences correctly retain the mutation; deduplicating event labels
would be wrong. If one true shared-mutant haplotype is withheld by a floor, its
truth event is a false negative and its allele remains missing. Cis/trans
assignments must continue to follow GT order, not majority/minority labels.

**Decision:** do not promote a support gate based on the phasing read-list count.
A future minimum-support policy would need separately validated assignment of
all retained primary records (for example haplotagging with explicit ambiguous
assignments), an independently specified adequacy criterion, and explicit
minority/no-call denominators. It still cannot recover input length classes
excluded before variant calling. Root is investigating that upstream boundary.

## Read-list integrity repair after independent audit

An independent review found that the experimental helper could adopt a phased
VCF alongside an empty tool read list. Twelve new failing regression cases
exposed this before the repair. The helper now validates the actual eight-column
read-list schema, known unique record IDs, single-input source ID, selected sample,
haplotype index, phase-set consistency, positive covered-variant count, and valid
first/last retained heterozygous positions. A phased result with no assignments
raises a visible error; an unresolved result with an empty list remains unresolved.
This repair adds no minimum support floor and makes no genotype changes.

After repair:27 deterministic unit tests and17 actual-tool integration tests
pass (44 total). Ruff, formatting and configured mypy pass. Five cached helper
reruns also retain their prior statuses:four phased and one not-needed. Their
fresh guard metadata is preserved in `production_cached/guard_verification/`.
Module SHA256 after this repair:
`947c4cb48613c9ca83835ec6ea52b0eb36650a6d5603cea02b68b0545e6cbedd`.
