I’ll inspect the implementation and frozen acceptance criteria, then check the benchmark evidence against each gate. I’ll use the review skills to distinguish verified findings from claims the available evidence cannot establish.


---

The supplied results already fail the ONT reconstruction gate: +1 exact diploid reconstruction against the required +5. The frozen specification also requires mutation sensitivity to stay at least at baseline, while the reported pooled sensitivity falls from 72/96 to 69/96. I’m checking the per-sample evidence and implementation, particularly the changing control denominator and whether haploid conversion is restricted to resolved alleles.


---

The control denominator changes because unresolved negatives are excluded from the reported specificity. The candidate has 11 confident negatives out of all 44 controls (25%), plus two false alarms and 31 unresolved negatives. That fails the frozen ≥65% confident-negative gate. I also found that missing alignment scores are counted as decisive read support, which can mistake suppressed secondary alignments for evidence of a second allele.


---

A stored development result contradicts the claim that homozygous false splits are eliminated: `ctrl_ident_80_80_hifi` is called as 80/75 and still raises a mutation alarm. The variant filter also caps every positive user QUAL threshold at 4 in haploid mode, so requesting QUAL 20 effectively applies QUAL 4. Read-only static checks pass; the full CI command is blocked because `uv` needs to create a lock file.


---

I found a benchmark input defect: every development FASTQ contains duplicate read names attached to different sequences. HiFi has 11,288 records but only 7,702 unique names; ONT has 13,400 records but only 8,287 unique names. The new dominance code merges scores by read name, so it can combine alignments from different reads. This undermines the support counts and requires an input-identity audit before the benchmark can support acceptance.


---

# Verdict: **REVISE AND RESUBMIT**

The candidate improves several development endpoints, particularly HiFi reconstruction and mutation-call precision. However, it fails multiple frozen acceptance gates, does not implement the required complete pairwise rescoring, and has a newly identified benchmark input defect that compromises read-level evidence.

**The protected 60 datasets should remain sealed.**

## Review scope

I inspected the working-tree implementation against `9535f7e`, including untracked modules, the frozen specification, Checkpoint 2 synthesis, both 140-sample evaluation reports, run configurations, development FASTQs, selected BAMs, and active VCFs.

I independently recomputed platform-specific results. All 140 baseline/candidate pairs have matching recorded input hashes; baseline reference hashes match the reference at `9535f7e`. Those are useful provenance checks, but they do not resolve the defects below.

---

## 1. Blocking findings

### 1.1 Every development FASTQ contains duplicate read names associated with different sequences

This is the most consequential finding.

| Platform | Datasets affected | FASTQ records | Unique read names |
|---|---:|---:|---:|
| HiFi | 70/70 | 11,288 | 7,702 |
| ONT | 70/70 | 13,400 | 8,287 |

For example, `ctrl_ident_80_80_hifi` contains **139 records but only 93 unique names**; 46 names identify different sequences.

The new score extractor stores:

```python
scores[qname][contig] = highest_alignment_score
```

Consequently, it can compare the score of **one physical read against c1 with another physical read against c2**, treating them as one observation. It also undercounts distinct reads. See [read_dominance.py:82](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/read_dominance.py:82).

This problem extends beyond dominance scoring: FASTQ conversion groups records by read name, and the existing extraction path uses `samtools fastq`. Its documented grouping and filtering behavior makes read-name uniqueness material to reconstruction. [Samtools documentation](https://www.htslib.org/doc/samtools-fasta.html).

**Required correction:** assign globally unique, neutral identifiers to simulated read records before mapping; verify sequence and quality preservation; regenerate input manifests and rerun **both** baseline and candidate. Merely adjusting the reported dominance counts cannot repair downstream artifacts already generated from colliding identifiers.

This does not establish that every observed improvement is spurious. It does mean the present benchmark cannot establish the intended read-level scientific validity.

### 1.2 Missing alignment scores are incorrectly treated as decisive evidence

In [read_dominance.py:131](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/read_dominance.py:131), a read observed only on c2 automatically increments `d2`.

That inference is invalid:

> Absence of a reported c1 alignment does not establish that `AS(c2) − AS(c1) > Δ`.

Minimap2 suppresses secondary alignments through score-ratio filtering and an output-count limit; the current mapping command does not request exhaustive candidate scores. The default secondary limit is five. [Minimap2 manual](https://lh3.github.io/minimap2/minimap2.html).

This was **explicitly identified and required to be fixed at Checkpoint 2**, which called for pairwise rescoring of spanning reads. The candidate instead reads the existing, censored ladder BAM.

A direct reproduction confirms that three records containing only `{"c2": 1}` are accepted as a second HiFi allele despite having no measured opposing scores.

**Required correction:** obtain comparable scores for the same eligible read against both candidates. Missing or incomparable scores must remain unresolved unless a justified bound establishes dominance.

### 1.3 Homozygous false splitting persists—and the minority rescue path causes it

The stored candidate result for `ctrl_ident_80_80_hifi` is:

- Allele 1: **80 repeats**
- Allele 2: **75 repeats**
- Mutation alarm: **present**

See its [alleles.json](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/tests/results/candidate_dev_140/ctrl_ident_80_80_hifi/alleles.json).

The actual cached-score comparison between `contig_71` and `contig_66` returns:

```text
d1 = 79
d2 = 3
ambiguous = 0
second allele accepted
```

**All three c2-supporting names lack a c1 score.** Thus the missing-score defect is exercised by a real false-positive development result.

The review request’s “100% rejected” statement is not true of final pipeline outputs. Rejecting an indel-valley split does not ensure specificity when a subsequent minority-rescue path can introduce another false split.

### 1.4 The validated candidate can differ from the emitted candidate

The rescue path evaluates the weighted cluster center, then `_with_support()` refines the selected contig again. The final output therefore need not be the contig whose dominance passed.

Additionally:

- Initial two-cluster discoveries bypass the new dominance validation.
- Empty score evidence can cause an indel split to be accepted.
- Candidate clusters can overlap existing cluster membership.
- The primary-support safeguard is not supplied by the production calls.

These are visible in [alleles.py:473](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/alleles.py:473).

**Required correction:** validate every final proposed pair, preserve the identities of validated candidates, and fail to an explicit unresolved state when required evidence is unavailable.

---

## 2. Numerical acceptance gates

The following results are independently recomputed from the stored evaluation reports. **“Pass” below means the reported artifact arithmetic passes; benchmark validity remains conditional on repairing read identities.**

| Frozen gate | HiFi | ONT | Assessment |
|---|---|---|---|
| Exact diploid sequence net improvement: ≥8 / ≥5 | 4→18; 14 wins, 0 losses | 9→10; 1 win, 0 losses | **HiFi passes; ONT fails** |
| HiFi exact reconstruction increase ≥12 percentage points | **+20.0 pp** | — | Pass |
| Exact count-pair wins > losses | 27→31; 4 wins, 0 losses | 15→17; 2 wins, 0 losses | Pass |
| No regression on previously exact count pairs | 0 losses | 0 losses | Pass |
| Mutation alarm sensitivity ≥baseline | 35/48→36/48 | **37/48→33/48** | **ONT fails** |
| Zero new control mutation alarms | 0 new | 0 new | Pass |
| Absolute control alarm cap ≤1 per platform | **2/22** | 0/22 | **HiFi fails** |
| Confident-negative control rate ≥65% | **5/22 = 22.7%** | **6/22 = 27.3%** | **Both fail** |
| Annotation FDR on mutation-positive samples decreases | 44.9%→25.6% | 39.5%→27.3% | Pass |
| Runtime increase ≤15%, specified repeated paired design | Evidence incomplete | Evidence incomplete | Not established |
| Peak process-tree RSS ≤8 GB | Not supplied | Not supplied | Not established |

Source reports: [baseline](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/tests/results/baseline_dev_140/evaluation_report.json), [candidate](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/tests/results/candidate_dev_140/evaluation_report.json).

### Exact annotation improvement is not demonstrated on the fixed truth denominator

Exact annotated-event recovery is:

- HiFi: **27/48→29/48**, an increase of **4.17 pp**.
- ONT: **26/48→24/48**, a decrease of **4.17 pp**.
- Pooled: **53/96→53/96**, no improvement.

Thus a ≥5 pp gate interpreted as exact recovery of true annotated events is not met. Annotation precision improves, but precision and recovery are different endpoints. If “annotation accuracy” was intended to mean another metric, its frozen numerator and denominator must be identified; it cannot be selected retrospectively.

Supported-event recovery improves from 32/96 to 43/96, but this is not interchangeable with raw exact annotation recovery or alarm sensitivity.

### The specificity denominator is conditional

The reported candidate specificity, **11/13**, comprises:

- 11 confident true negatives;
- 2 false-positive controls;
- **31 unresolved negative controls excluded from that denominator**.

Across all 44 controls, the confident-negative rate is **11/44 = 25%**.

The remaining false alarms are:

- `ctrl_ident_80_80_hifi`
- `ctrl_lowcov_60_80_hifi`

Both were already false positives in baseline. Therefore **zero new alarms passes**, while **zero absolute alarms fails**, as does the frozen HiFi cap of one.

### “Zero regressions” needs a precise scope

There are zero losses in previously exact **whole-sample sequence reconstructions** and exact **count pairs**.

There are nevertheless mutation regressions. For example:

- `mut_eq_inscccc_60_60_h1_hifi`: exact event TP 1→0.
- `mut_asym_ins16bp_20_90_h2_ont`: TP 1→0.
- `mut_gap1_insg_50_51_h1_ont`: TP 1→0.
- `mut_gap3_dupc_60_63_h1_ont`: TP 1→0.

Missing allele observations also increase from **39 to 60**. Some reduction in output cardinality correctly removes false splits; the aggregate cannot establish that every removed observation was inappropriate.

### Statistical interpretation

For the primary binary endpoint, nominal two-sided exact McNemar tests give:

- HiFi, 14 wins/0 losses: **p = 0.000122**.
- ONT, 1 win/0 losses: **p = 1.0**.

These are exploratory development results, subject to candidate selection and the input defect. They do not provide confirmatory significance after tuning.

Also, the report’s `min == max` bounds concern alternative sequence assignments. **They are not confidence intervals.**

---

## 3. Scientific validity of the mechanisms

### 3.1 Comparing secondary candidates with the primary peak

**Conceptually appropriate; insufficient as implemented.**

Using the best-supported primary contig as the comparator is preferable to validating a secondary candidate against an arbitrarily weakened alternative. It does not inherently leak simulation truth.

However, it does not independently validate:

- whether the primary peak is correct;
- whether both scores exist for each read;
- whether reads span the informative interval;
- whether identifiers represent distinct reads;
- whether final read partitions are pure.

The primary peak itself is chosen from conditional alignment summaries, with different reads contributing to different contigs. Calling it “validated” overstates the evidence.

Candidate discovery and validation may legitimately use the same reads, but searching multiple secondary candidates requires empirical false-positive calibration. Changing the comparator alone provides no specificity guarantee.

### 3.2 The 1% ratio and absolute read floor

**Reasonable candidate heuristics to investigate; not established biological or statistical thresholds.**

The implemented requirement is approximately:

\[
d_2 \ge \max\left(m,\left\lceil0.01(d_1+d_2)\right\rceil\right)
\]

with **m=3 for HiFi and m=4 for ONT**, contrary to the request’s platform-independent “≥3” wording.

Three concerns matter:

1. **The denominator excludes ambiguous reads.**  
   I reproduced acceptance with three c2-dominant reads and 1,000 ambiguous reads: the reported dominance ratio is 100%, although only 0.3% of scored reads are decisive for c2.

2. **Three read names do not establish three independent molecules.**  
   This is especially problematic here because identifiers collide. PCR artifacts, chimeras and systematic errors also invalidate a simple independent-read interpretation.

3. **There is no universal false-positive control.**  
   The rule depends on depth, candidate-search breadth, read quality, completeness and error correlation. Extreme PCR asymmetry motivates permissive discovery, but does not itself justify accepting an allele.

The margins 43/42 are derived from an assumed 60-base gap penalty. That is a mechanistic heuristic, **not empirical calibration**. Explicit minimap2 preset overrides can also change score interpretation without changing these platform-keyed thresholds.

Required evidence includes tuning/selection separation, threshold sensitivity, eligible-read counts, ambiguity fractions, and false-positive performance under null homozygous data.

### 3.3 Haploid majority re-genotyping

**Conditionally justified for a pure, resolved haploid partition; not justified merely because a BAM maps to one reference contig.**

Under a simple biallelic model with independent symmetric errors, majority voting is defensible. Encoding the chosen allele as `1/1` can serve as an internal consensus convention. `bcftools consensus` uses the selected genotype and does not perform this inference itself. [Bcftools consensus documentation](https://samtools.github.io/bcftools/howtos/consensus-sequence.html).

The implementation’s central assumption is unverified:

- Dominant read identities are not passed into calling.
- Calling still extracts reads through cluster-contig membership.
- `samtools fastq` excludes secondary alignments by default, so reads counted as supporting discovery can disappear before consensus.
- Distinct inferred lengths do not establish clean haplotype separation.

The explicit same-length route remains diploid, which is a positive feature. But an incorrectly split equal-length sample bypasses that protection.

Additional defects in [vcf.py:82](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/vcf.py:82):

- **User QUAL thresholds are silently weakened:** `min(min_qual, 4.0)` turns QUAL 20 into QUAL 4.
- **Multiallelic AF is unsupported:** `float("0.5,0.4375")` fails and silently retains the original genotype.
- I found **15 multiallelic heterozygous records in active distinct-length candidate VCFs**. The purported elimination of ambiguity is incomplete.
- The code rewrites all parseable genotypes, including homozygous and potentially missing calls.
- AF=0.5 ties are forced ALT.
- Genotype confidence fields are retained despite the changed genotype interpretation.

A concrete active record in `ctrl_asym_25_100_ont` remains:

```text
GT=1/2   AF=0.5,0.4375
```

Required correction: explicit partition-purity evidence, allele-aware handling of multiallelic records, documented tie/missing-data policies, preserved original calls, and an explicit configurable consensus threshold that respects user overrides.

---

## 4. Reference architecture and claimed causal evidence

The proximal slicing change is internally coherent: **all 150 bundled contigs exactly match current `build_contig()` output**.

Nevertheless, two claims need correction.

### “All 44 boundary failures eliminated” is unsupported

The baseline failure atlas contains 44 entries labeled `variant_calling`, but:

- Only **17** belong to its `boundary_mutation` category.
- The category includes missed or misannotated events; it does not demonstrate Clair3 execution failures.
- **Six** of those 44 samples still have event false negatives in the candidate report.

The current diagnosis script also tests sequence discordance before mutation failure. A change in stage labels cannot demonstrate elimination of the underlying errors.

A causal attribution to the flank change requires a matched flank-only ablation and per-site evidence, separating it from majority conversion, QUAL changes and mutation-template changes.

### Coordinates and reference compatibility need resolution

The request cites `chr1:155160475–155160974`, while the bundled dictionary identifies the VNTR as `chr1:155188487–155192239`. These must be reconciled with an explicit assembly and coordinate convention.

The compatibility guard also omits `proximal_flank`: setting it false can retain the bundled proximal reference while selecting distal trimming anchors. See [pipeline.py:63](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/pipeline.py:63).

Exact anchor matching can fall back to fixed trimming when anchors are altered or displaced. Report and test those cases; reference-filled flanking sequence must not be mistaken for observed read support.

---

## 5. Architecture, challenging strata and overfitting

### Modular separation is incomplete

`discover_length_candidates()` and `validate_candidates_with_dominance()` are tested but are **not used by production `detect_alleles()`**. Production imports only the evidence-reporting helper from that module.

Consequently, high coverage of `length_candidates.py` does not validate the actual rescue implementation.

The new margins, read floors, ratio, AF threshold and QUAL cap are also not exposed through typed immutable runtime settings. This fails the frozen configuration contract.

### Challenging-stratum gates remain problematic

- **Identical controls:** single-output/no-alarm behavior occurs in 6/7 HiFi and 7/7 ONT cases. If “correct single sequence” requires exact sequence, ONT achieves only **5/7**. Neither interpretation supports the blanket 100% claim.
- **HiFi gap 1:** all six development gap-1 datasets have `all_counts_exact=0`. The required ≥20 usable reads per allele eligibility table is missing, so the conditional gate cannot be formally adjudicated. The observed results provide no evidence of success.
- **25/140 asymmetry:** some HiFi recoveries occur, but both ordinary ONT cases remain count-inexact. The ≥5 usable minority-read denominator is unreported and currently compromised by identifier collisions.
- The low-support minority no-call contract requires explicit per-case validation; aggregate reconstruction metrics do not establish it.

### Benchmark protection and provenance need strengthening

I found no explicit truth-name or simulation-label branching in the new runtime modules. However:

- The evaluation driver supplies `raw_reads_path` and descriptive output names containing lengths and mutations, bypassing the neutral-token execution route.
- It reads the sealed ledger before filtering by split.
- The generator prints the ledger hash **after generation**, rather than demonstrating the specified pre-generation commitment.
- The stated seed non-overlap requirement is ≥1,000, while implemented seeds advance by **10**.
- The claimed inner 50/20 tuning/selection evidence is absent from the supplied comparison.

These are protocol weaknesses, not proof of deliberate leakage.

There is also an unreported scientific change in [config.py:71](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/config.py:71): `delete_insert` now retains both boundary bases, changing the bundled `delinsAT` template from a net −1-base to a net +1-base edit. It may be a valid correction, but needs independent sequence-level justification, focused tests and explicit evaluator provenance. Agreement with the simulator alone is insufficient.

---

## 6. Required resubmission

Before opening final validation:

1. **Repair read identifiers and audit simulation invariants**, then rerun both comparators with fresh artifacts and manifests.
2. **Implement complete pairwise scoring of eligible spanning reads**, preserving physical-read identity and treating missing scores as unresolved.
3. **Validate the final emitted candidate pair on every discovery path**, and carry validated read assignments into calling.
4. **Make majority consensus partition-aware and multiallelic-aware**, preserve original evidence, and respect configured quality thresholds.
5. **Provide controlled ablations** for flank, dominance, majority conversion, QUAL and mutation-template changes.
6. **Reproduce every frozen gate by platform and eligible stratum**, including mutation regressions, unresolved controls, repeated runtime measurements and process-tree RSS.
7. **Freeze code, resources, evaluator and parameters** only after development acceptance and protection-protocol corrections.

## Verification performed

Read-only verification passed:

- Ruff lint and formatting: **107 files formatted**.
- Mypy: **no issues in 53 source files**.
- File-size gate: **123 files checked**.
- Actionlint.
- Four existing pure dominance tests.
- Independent report arithmetic, identifier audit and selected artifact checks.

`make quality` could not run through `uv` because the read-only environment prevents lock-file creation; the underlying static checks were run directly. I did **not** independently reproduce the claimed 707-test coverage run, full integration suite or packaging checks. No files were changed, and protected validation datasets were not inspected or executed.

**Final recommendation: REVISE AND RESUBMIT.** The HiFi improvement warrants continued development, but failed ONT and control gates, invalid read identity, incomplete rescoring and unsafe consensus assumptions prevent proceeding to confirmatory validation.