# Equal-length false split: bounded development investigation

2026-09-14. Production code and tests were not changed by this investigation.
Stopped after the user's instruction to finish the current work and release.
No final/reserved seeds, mapping, new calling, or parameter promotion occurred.

## Decision

The two experimental span-mixture models are **not promoted**. Both still split
true equal 60/60 into false 59/60. They improve exploratory distinct-length sets
on other exposed inputs but do not fix the selected failure, establish molecular
haplotypes, or validate the required downstream mutation gate. The release must
retain the equal-length failure disclosure. No test was weakened.

## Cached BAM reproduction and cause

Calling the current `detect_alleles` on cached
`development/hifi/sample_homozygous_60_60/mapping.bam` reproduces 60/63.
The BAM has 169 primary records and 1,014 alignment records. The selected
60-repeat cluster has 723 alignment records/164 primary records; the false
63-repeat cluster has 276 alignment records/two primary records. These counts
are alignment observations, not independent PCR molecules.

| Selected reference | Alignment records | Primary records on this reference | Mean AS | Mean indel bp | Mean clipped bp | Covers both reference VNTR boundaries |
|---|---:|---:|---:|---:|---:|---:|
| contig_51 (60 total) | 167 | 0 | 3170.07 | 25.11 | 5.16 | 158 |
| contig_54 (63 total) | 111 | 0 | 3093.49 | 34.74 | 51.89 | 1 |

The second low-indel valley mostly comprises alignments that do not traverse the
longer reference VNTR. Indel-only means therefore mix unlike alignment extents;
low indel burden does not prove a complete longer allele. Primary-only selection
is also insufficient: the correct reference has zero primary alignments, while
contig_49 has 125. A primary designation is not a per-read maximum-AS assertion.

For the per-read comparison, QNAME plus full query length was used **only** where
exactly one primary record exists. This identifies 161 groups. Four colliding
name/length groups contain eight primary records and are explicitly ambiguous.
Among 103 unambiguously paired records with alignments on both selected
references, every AS difference `AS(contig_54)-AS(contig_51)` is negative,
ranging from -150 to -85. No paired record prefers the false longer candidate.

This motivates an unimplemented candidate-dominance diagnostic using comparable
per-record fit scores. It is not a validated likelihood: minimap2 censors the
reported secondary alternatives; absent counterfactual alignments cannot be
assigned invented scores; ambiguous identities cannot be arbitrarily joined;
and AS differences require calibration before probabilistic interpretation.
No support floor or clipping threshold was introduced.

## Predeclared bounded span models

The existing `read_evidence/anchors.json` inventory supplies all 77 exposed jobs:
44 original HiFi, three ONT, six previously exposed controls, and 24 dependent
read perturbations. Models use unrounded observed span midpoints in base pairs
at exact/one-edit/two-edit anchor settings. Truth is consulted only afterward
for scoring. No missing input is dropped from the denominator.

Two alternatives were executed: Gaussian and fixed-df4 Student-t mixtures,
each comparing one versus two components with one shared variance. BIC parameter
counts are two and four respectively. Four fixed quantile starts (1/99, 10/90,
25/75, 40/60 percentiles) allow minority components to be considered. No outlier
removal, minimum-read filter, mode suppression, or jitter was used. Shared variance
prevents a separate component variance collapsing around a singleton, but is an
imperfect assumption for length-dependent errors. A numerical variance floor of
1e-9 was never reached; all retained fits converged within the 200-iteration cap.

These are exploratory mixture scores, not calibrated genotype confidence.
Mixture nulls are nonregular, sequencing errors need not follow either assumed
distribution, and midpoint spans do not model boundary uncertainty or partial
read ascertainment. Rounded fitted means are only a distinct-length diagnostic;
a one-component result cannot establish two identical molecular alleles.

At exact and one-edit anchor settings, both models yield:

| Cohort | Current exact distinct-length sets | Model exact distinct-length sets | Previously exact sets lost |
|---|---:|---:|---:|
| HiFi | 40/44 | 43/44 | 0 |
| ONT | 0/3 | 3/3 | 0 |
| Exposed controls | 2/6 | 5/6 | 0 |
| Perturbations | 9/24 | 23/24 | 0 |

At two edits the Gaussian model falls to 1/3 ONT; Student-t remains 3/3.
All other cohort totals are unchanged. Real 60/61, 60/62, and 60/63 controls
remain distinct. The 25/140 minority component is retained. **Both equal-length
60/60 controls and the renamed equal-length perturbation still fail**, with
59/60 selected at every tested setting.

For the original equal-length sample at one-edit anchors, the Gaussian component
means are 3557.68 and 3604.67 bp, with BIC improvement 94.37 for two components;
Student-t means are 3556.27 and 3605.09 bp, with improvement 35.20. Thus replacing
indel valleys with BIC still treats the short error tail as an extra length.
A larger arbitrary separation threshold would threaten the real 60/61 case and
was not tried. The 24/28 downstream mutation acceptance requirement was not
re-evaluated because no pipeline candidate was implemented or called.

## Primary-method context

The prior GeneFoundry/PubTator METHODS retrieval is documented in
[the research note](genefoundry-pubtator-research.md) and
[the length-model check](length-model-research.md). HMMSTR explicitly describes
GMM overcalling narrow nearly homozygous distributions, consistent with this
failed prototype. Straglr's mixture approach and NanoRepeat's rounding/jitter
caveat motivate controls, not automatic transfer to MUC1.

- Van Deynze K et al. Enhanced detection and genotyping of disease-associated tandem repeats using HMMSTR and targeted long-read sequencing. Nucleic Acids Res. 2024. PMID:39676678. doi:10.1093/nar/gkae1202. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11754662/).
- Chiu R et al. Straglr: discovering and genotyping tandem repeat expansions using whole genome long-read sequences. Genome Biol. 2021. PMID:34389037. doi:10.1186/s13059-021-02447-3. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8361843/).
- Fang L et al. Haplotyping SNPs for allele-specific gene editing of the expanded huntingtin allele using long-read sequencing. HGG Adv. 2023. PMID:36262216. doi:10.1016/j.xhgg.2022.100146. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9574884/).

## Artifacts and executed commands

Ignored directory: `tests/results/production_validation_20260914/length_failure_fix/`.
`equal_bam.py` / `equal_bam.json` retain per-contig metrics, paired records, and
identity ambiguities. `equal_models.py` / `equal_models.json` retain the explicit
77-input results, model settings, BICs, fitted parameters, convergence flags,
source-evidence hashes, and all offline scores. No new dependencies were installed;
the attempted NumPy availability check found it absent, so the prototype uses
Python's standard library.

```sh
uv run --locked --all-extras python tests/results/production_validation_20260914/length_failure_fix/equal_models.py --pv tests/results/production_validation_20260914 --output tests/results/production_validation_20260914/length_failure_fix/equal_models.json
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH uv run --locked --all-extras python tests/results/production_validation_20260914/length_failure_fix/equal_bam.py --pv tests/results/production_validation_20260914 --output tests/results/production_validation_20260914/length_failure_fix/equal_bam.json
```

Both scripts refuse existing output paths. These are unfinished research
prototypes, not maintained production implementations or a release feature.
