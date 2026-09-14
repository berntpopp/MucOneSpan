# Deep MucOneUp validation of MucOneSpan 0.10.0

Date: 2026-09-14. Baseline commit: `d8390b3`. Production algorithms were held
fixed throughout the baseline. Experimental replacements ran only in separate
processes/output directories. No changes have been committed or released.

## What the fresh baseline establishes

The current pipeline detects many known mutations and usually reports the right
repeat counts, but it does **not** reliably reconstruct both full VNTR alleles.
On the 44 existing HiFi samples, both complete nucleotide sequences were exact
in **10/44 samples**, despite exact repeat-count pairs in **40/44**. Detection
of any mutation and exact mutation identification are materially different.

| Endpoint | Existing HiFi challenge set | Newly generated ONT challenge set |
| --- | ---: | ---: |
| Pipeline completion | 44/44 | 3/3 |
| Correct mutation parent, name, repeat index and assigned haplotype | 24/28 | 2/2 |
| Any mutation reported in a mutant sample | 25/28 | 2/2 |
| Normal samples with no mutation reported | 15/16 | 1/1 |
| Exact reported repeat count, individual alleles | 82/88 | 3/6 |
| Exact pair of reported repeat counts | 40/44 | 0/3 |
| Both reported counts within ±2 repeats | 40/44 | 3/3 |
| Exact complete nucleotide sequence and repeat structure, individual alleles | 39/88 | 3/6 |
| Both complete sequences exact | 10/44 | 0/3 |
| Extra mutation records, including extras on positive samples | 14 | 0 |
| Extra frameshifts marked VCF-supported | 7 | 0 |

Exact mutation sensitivity is **85.7%**, with a descriptive Wilson 95% interval
of **68.5–94.3%**; normal-sample specificity is **93.8%**, interval
**71.7–98.9%**. The historical inclusive detection rate, 25/28 = 89.3%, counts a
wrong/partial mutation identification as detected. Exact event precision is
24/(24+14) = **63.2%** across all reported mutation records, or 24/(24+7) =
**77.4%** among VCF-supported frameshifts. Sample specificity alone conceals
additional wrong mutations on true-positive samples.

These are descriptive challenge-set statistics. Most positives are dupC,
several seeds reuse related designs, and the later subsamples overlap. They
are not independent population performance estimates. ONT has only two
positives and one negative; its nominal 100% detection/specificity has very wide
Wilson intervals (34.2–100% and 20.7–100%, respectively). No confidence interval
treating two alleles of the same sample as independent is claimed.

The exact HiFi mutation failures are dupC at repeat 45 in 100/120, insCCCC at
repeat 25 in 60/80, dupC at repeat 25 in the same-length 60/60 pair, and dupC at
repeat 50 in 120/140. The normal-sample false positive is 40/44, with five
reported mutations, two marked VCF-supported. Complete per-allele results,
including errors and false extra calls, are in the accuracy JSON artifacts.

## Truth and experimental validity

- All 88 original haplotype FASTAs exactly reconstruct from the dictionary,
  stored mutant unit and explicit flanks. All 28 mutated units match the
  simulator mutation definitions. With new ONT truth, direct classification is
  exact for **94/94 haplotypes**; **132/132** supported mutation-template/context
  probes also pass. This isolates much of the reconstruction loss upstream of
  classification on ideal sequences.
- Repeat counts include five pre-repeat and four after-repeat units.
  `contig_N` contains N canonical units, so its total count is N+9. The benchmark
  uses stored structures/stats, not sample filenames. The historical benchmark
  incorrectly describes all five `sample_bench_*` cases as 60/80.
- Stored mutant `repeat_lengths` arrays incorrectly retain 60 for mutated
  units; sequence truth is authoritative. All original metadata coverage
  fields say 30 despite generation commands requesting 200 or 50. Current
  simulator help defines amplicon coverage as initial template molecules,
  not achieved read depth. Actual usable reads and allele balance must be
  measured; the asymmetric sample contains very few long molecules.
- Original simulator version is 0.44.2; new ONT generation uses current local
  MucOneUp 0.44.5 and QSHMM-ONT-HQ. These platforms and simulator versions are
  reported separately. ONT simulation used mean accuracy 0.95; this does not
  establish performance on every chemistry/basecaller or Q20+ dataset.
- Simulator read names collide across distinct haplotypes. Record counts are
  used here rather than treating unique QNAMEs as unique molecules. A direct
  check found no FASTQ record loss in the tested baseline conversions; name
  collisions are nevertheless unsuitable input for name-based phasing.
- Pipeline truth matching uses a one-to-one permutation minimizing total exact
  global sequence edit distance. All original HiFi/ONT assignments had a unique
  minimum. This is a scoring procedure, not proof of read-backed phasing.
- Inputs, truth, dictionary, reference and model fingerprints, exact CLI
  arguments, versions, logs and per-stage times are retained in ignored results.

## Causes established by targeted experiments

### Secondary alignments distort coverage, but removing them damages length inference

The false second 60/60 allele is reported with 276 alignment records but reaches
Clair3 with only two primary reads. The 25/140 sample instead reports 25/28 and
feeds the spurious 28-repeat branch one primary read. Secondary mappings are
useful alternative-reference fit evidence; they are not independent support.

An ablation on all **44 fresh mapping BAMs** confirms that simply counting only
primary alignments changes exact pair recovery from **40/44 to 0/44**, commonly
shifting both lengths down by 1–2 units. Keeping primary-only refinement does
not fix that. Separate candidate scoring from molecule support rather than
changing the histogram filter globally. Current per-contig coverage thresholding
also removes diffuse allele evidence before pooling, and flat indel plateaus
can trigger an arbitrary valley split.

### Default consensus creates ambiguity; explicit ALT is an improvement with tradeoffs

The installed bcftools 1.17 selects genotype-based IUPAC consensus by default
for these sample-containing VCFs. Source inspection and a two-SNP phased
reproduction establish this behavior, contradicting the pipeline comment that
all ALT alleles necessarily become an unambiguous mutant haplotype.

Holding raw Clair3 calls and QUAL=5 fixed, explicit `-H A` changes exact HiFi
sequence recovery from **39/88 to 60/88**, and complete diploid reconstruction
from **10/44 to 24/44**. It removes 1,189 ambiguous bases. Exact mutation TP
stays **24/28**, while VCF-supported extra calls rise **7 to 13**. ONT exact
reconstruction remains 3/6 despite removing 107 ambiguous bases. This is useful
evidence for specifying consensus behavior, not justification for a blind
default switch.

For equal-length alleles, the pipeline splits homozygous-ALT sites into one
VCF and all sites into another. Opposite-phase heterozygous variants cannot be
reconstructed correctly by that rule, even with explicit ALT selection. A
real bcftools probe with two trans SNPs produces reference and IUPAC mosaics
instead of the two known haplotypes. Phase or report unresolved structure.

### The long dupC failures have different mechanisms

The 100/120 raw VCF contains the truth-consistent insertion at contig_91:3192,
QUAL 3.22, PASS. QUAL>=5 removes it. The 120/140 case has many SNP records but
no raw indel. Both contradict the historical blanket explanation of zero
Clair3 variants; only the former is demonstrated to be threshold loss. The
cached-VCF threshold sweep is documented separately below. Detecting a motif
somewhere in a repetitive read is insufficient evidence to rescue the latter.

All ten cached-call conditions (five QUAL cutoffs × two consensus modes) ran
on all 47 samples, **470 reconstructions** in total. Selected HiFi comparisons:

| QUAL / consensus | Exact sequences / 88 | Exact diploid / 44 | Exact mutations / 28 | Supported extra calls | Normal FP / 16 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 5 / default baseline | 39 | 10 | 24 | 7 | 1 |
| 5 / explicit ALT | 60 | 24 | 24 | 13 | 1 |
| 3 / explicit ALT | 66 | 27 | 26 | 13 | 1 |
| 20 / explicit ALT | 14 | 1 | 14 | 0 | 0 |

QUAL3/ALT recovers the 100/120 dupC and 60/80 insCCCC (QUAL 4.47), but adds a
supported false positive in the sole ONT normal sample. QUAL20 sacrifices half
the true HiFi mutations and most exact structures. The best result on this
challenge set is not a validated platform-independent operating threshold.

### Several confidence and correctness contracts need repair

- A 59-base repeat with two deleted bases and one insertion is marked
  non-frameshift because the classifier sums absolute indel lengths instead
  of signed net change.
- The advertised backward rescue requires >30 unconsumed bases, while normal
  forward processing consumes until fewer than 30 remain. A 40-base insertion
  produces two wrong mutation calls without reaching the intended rescue.
- Five perfect repeats followed by 29 Ns still receive confidence 1.0 and
  100% exact matches; unclassified sequence is omitted from the score.
- `vcf_support` accepts any nearby variant position, including a SNP in the
  following repeat, without checking REF/ALT or the named mutation. It is
  positional proximity, not independent validation of the indel.
- A consensus with 286 nucleotide edits and 59 repeat-structure edits receives
  average confidence **0.9828**. The score is not a calibrated probability of
  exact reconstruction.
- VCF query failures can become empty variant lists, and absent heterozygous
  calls can become a homozygous designation. Error, insufficient evidence,
  equal length and identical sequence need separate states.

## Performance measured

The 44 original HiFi pipelines took **453.3 seconds** in summed wall time,
median **6.52 s**, maximum **36.36 s**, using four requested threads per run.
Other audit jobs ran concurrently, so these are local measurements rather than
machine-independent throughput guarantees.

| Stage | Summed time | Share of wall time |
| --- | ---: | ---: |
| Mapping | 47.88 s | 10.6% |
| Allele detection | 0.61 s | 0.1% |
| Calling, including extraction/remapping/filtering | 232.84 s | 51.4% |
| Consensus construction | 1.32 s | 0.3% |
| Classification | 169.11 s | 37.3% |

One difficult consensus makes 80,240 scalar edit-distance calls; those consume
96.5% of its profiled runtime. Three alternating unprofiled comparisons with an
exact Python bit-vector implementation reduce its median from **19.350 s to
2.397 s**, **8.07× faster**, with identical complete return dictionaries. An
exact 120-repeat control remains about 0.58 ms. The full default-QUAL5 replay
also reproduces every baseline consensus and validated classification on all
47 samples. This is an experimentally supported optimization opportunity with
preserved scoring semantics, not an 8× whole-pipeline speed claim.

Per-allele remapping still buffers full FASTQ/SAM strings, repeats scans/indexing,
and can exceed the nominal thread budget across concurrent tasks. Streaming,
reference reuse and coordinated thread allocation are plausible further
optimizations; their speed/RSS effects were not measured in this study.

## Additional experiments and rejected shortcuts

**Depth:** eighteen deterministic, record-based subsets of the baseline dupC
and normal 60/80 inputs use 10, 20 or 40 reads, each with three seeds. All six
10-read runs produce explicit coverage errors. At 20 reads, dupC is found in
3/3 subsets but the complete two-allele structure is not recovered; one run
omits the second haplotype. At 40 reads, dupC is exact in 2/3, and 2/3 normal
subsets produce false mutations. Counts alone do not justify a confident
complete genotype at this depth. These overlapping subsets are not independent
sensitivity/specificity replicates.

**Read names:** six reruns change only QNAMEs to unique record IDs; the
sequence/quality multisets are verified identical. The baseline dupC 60/80,
long dupC 120/140, and normal 60/80 classifications change, while the original same-length
and long-mutation failures remain. The renamed normal sample acquires two
VCF-supported false calls. Read identity affects the mapping/calling path;
renaming is necessary input hygiene for phasing, not an accuracy fix by itself.

**Direct whole-template scanning:** scan all 44 known mutated templates on both
strands, counting each mutation name at most once per read. At >=3 supporting
reads it detects the expected name in 28/28 HiFi mutants, but falsely flags
16/16 normal samples and reports 96 extra names on positives. Increasing
support to ten reads still flags every normal. A deletion template occurs as
an internal substring of an ordinary repeat; sequencing errors also generate
mutation-like motifs elsewhere in long arrays. Reject this unlocalized rescue.

**Read-anchor length prototype:** a fixed rule selects supported span modes
between terminal-repeat anchors, without using truth to choose peaks. On the
original inputs it recovers all distinct length sets in 44/44 HiFi and 3/3 ONT
at >=3 exact-bin reads. Raising support to five loses the asymmetric minor
allele. It rejects 81.3% of ONT reads because exact anchors are absent. Its
±1-bin suppression cannot resolve a one-repeat separation, and a single mode
does not prove homozygosity.

**Held-out seeds:** six new HiFi samples use seeds 9101–9106, current MucOneUp
0.44.5 and 200 requested template molecules. Selection rules and anchor support
threshold three were fixed before these reads were generated. These samples
test gaps one/two/three, equal lengths, a mutation on the longer second
haplotype and a mutation near the terminal canonical repeat.

| Truth case | Current pipeline reported count pair | Exact mutation result | Fixed anchor distinct-length result |
| --- | --- | --- | --- |
| Normal 60/61 | 61/64 | No false calls | 60 only; misses 61 |
| Normal 60/62 | 61/61 | No false calls | 60/62 |
| dupC 60/63, haplotype 1 repeat 25 | 61/66 | Missed | 60/63 |
| Normal 60/60 | 60/63 | No false calls | Single length 60 |
| dupC 60/80, haplotype 2 repeat 25 | 60/80 | Exact | 60/80 |
| insCCCC 60/80, haplotype 1 repeat 55 | 60/80 | Exact | 60/80 |

All six full pipelines complete. Only **1/12 full allele sequences** and **2/6
count pairs** are exact; mutations are exact in **2/3**, with **0/3 normal
false-positive samples**. The fixed anchor prototype recovers **5/6 distinct
length sets**, with the predicted one-repeat-gap failure and no extra modes.
It is a length-evidence experiment, not a phased VNTR reconstruction. Direct
classification of the twelve new truth haplotypes is exact, bringing the
truth-only sequence checks to **106/106** across 53 full-data simulations.

Overall, **77 full pipeline attempts** were executed: 53 full-data simulations
and 24 depth/name perturbations. **71 completed; six 10-read subsets failed
with explicit coverage errors.** Cached consensus experiments are counted
separately and do not represent another 470 independent sequencing samples.

## Recommended implementation order

1. Correct benchmark endpoints, truth validation and explicit no-call states.
   Keep exact mutation identity/position, false extras and full structure in
   acceptance criteria; repair signed frameshift accounting and propagate
   external-tool failures with focused regression cases.
2. Replace scalar edit-distance scoring with the experimentally equivalent
   bit-vector implementation; retain traceback/tie behavior and add adversarial
   equivalence checks. This offers a measured speed gain without changing
   scientific thresholds or adding a dependency.
3. Specify consensus sample/haplotype behavior and investigate false extras
   before adopting explicit ALT. Preserve phase for equal-length alleles;
   expose ambiguity and unclassified sequence in reports and confidence.
4. Redesign length support around spanning molecules and robust anchor/sequence
   evidence while retaining secondary candidate-fit information. Validate on
   independent seeds, one-/two-repeat gaps, minority alleles, ONT anchor errors
   and partial reads. Do not silently force two lengths or homozygosity.
5. Rescue low-quality mutations only with independent, localized allele
   evidence. Use per-locus read support, strand/balance, variant identity and
   correct coordinate projection; evaluate extra calls and no-calls alongside
   rescued positives. A global QUAL reduction or motif-anywhere rule is not
   sufficient.
6. Then measure streaming, reference caching and thread budgets on larger
   molecule counts, with output-equivalence and peak-memory checks.

## Evidence and reproduction

Main ignored artifact root: `tests/results/deep_validation_20260914/`.

- `run_matrix.py`: unmodified CLI execution, stage timing, input hash, log and
  exit status per sample; refuses existing sample output directories.
- `provenance.json`, `hifi/measurements.json`, `ont/measurements.json`:
  versions, model/truth fingerprints and exact run arguments.
- `{hifi,ont,perturbation}_accuracy.json`, `statistics.json`: exact scoring,
  matched/missing alleles, false extras and descriptive statistics.
- `make_perturbations.py`, `perturbations.json`, `rename_comparison.json`:
  fixed subsample seeds and name-only checks.
- `direct_templates.py`, `direct_template_results.json`: rejected motif screen.
- `mapping_audit/`: fresh/historical primary-only ablations, anchor prototype,
  held-out generation commands and truth-independent mode selection.
- `consensus_experiments.py`, `consensus_experiments/`: raw-call-fixed QUAL and
  explicit-ALT experiments, plus baseline replay equivalence assertions.
- `results/classification-audit/`: truth probes, exact evaluator, adversarial
  cases, scalar/bit-vector equivalence and paired profiling.

Detailed companion reports:
[truth audit](2026-09-14-truth-audit.md),
[mapping and allele audit](2026-09-14-mapping-audit.md),
[classification audit](2026-09-14-classification-audit.md), and
[consensus/threshold experiments](2026-09-14-consensus-experiments.md).

External behavior was cross-checked against the primary
[minimap2 documentation](https://github.com/lh3/minimap2),
[Clair3 documentation](https://github.com/HKU-BAL/Clair3), and
[bcftools 1.17 consensus source](https://github.com/samtools/bcftools/blob/1.17/consensus.c).
Installed versions, rather than current upstream defaults, define this run:
minimap2 2.28-r1209, samtools 1.15.1, bcftools 1.17, Clair3 1.0.10, Python 3.12.9.

Checks actually executed: `make dev`; `make ci-check` (**239 passed, 86.11%
branch-aware coverage**, all static gates passed); tool-enabled `make test-int`
(**10 passed**, eight end-to-end cases deselected); explicit end-to-end selection
(**6 passed, 2 expected failures**, ten tool-integration cases deselected).
The two expected failures are same-length 60/60 and asymmetric 25/140 count
recovery. There were no missing-tool/model/data skips in these executions.
`make docs-check` also passed after adding the study link to the benchmarking
guide. The final diff contains reports/documentation only; packaging, security
and additional clinical-data validation were not part of these changes.
