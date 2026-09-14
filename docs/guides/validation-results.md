# Validation evidence for 0.11.0

This release improves exact scoring, evidence handling, evaluation and runtime
configuration. It does **not** establish a general improvement in complete diploid
MUC1 reconstruction. Two known length-detection failures remain. The planned new
32-simulation final panel was not generated: the release was scoped to the verified
work and a reproducible workflow for future experiments.

## What was verified

| Change | Evidence | Interpretation |
| --- | --- | --- |
| Exact bit-vector edit distance | At least 10,000 independent scalar-oracle comparisons; 141 cached complete classification outputs identical | Scientifically equivalent scoring optimization |
| Classification performance | Three alternating paired difficult-panel repetitions: median 17.612304 s to 2.214009 s | 7.955× local classification speedup, not a pipeline speedup |
| Default configuration integration | Five actual before/after cases, ten successful full pipeline runs across HiFi and ONT | Sequences, structures, counts and mutation records unchanged; provenance added |
| Corrected boundary anchors | 137 available cached full consensuses unchanged; 17 missing allele slots accounted for | Fixes anchor metadata; a synthetic flank insertion demonstrates the intended trimming correction |
| Signed indels and variant support | Deterministic mixed-indel tests and real bcftools normalization/consensus fixtures | Frame uses net inserted minus deleted bases; support requires exact projected identity |
| Failure and ambiguity accounting | Adversarial evaluation tests for wrong identity/position/parent, tied assignments, extra/missing alleles, stale files and failed commands | Failed or ambiguous samples cannot become successful negatives |

The configuration comparison is a functional regression check. Its index
warm/cold conditions were not controlled for a performance claim. Whole-pipeline
speed improvement has not been established. Confidence values remain heuristic
scores, not calibrated probabilities. VCF concordance is agreement with the VCF
used for consensus, not independent read validation.

## Development accuracy

All 44 original HiFi samples, three ONT samples, six subsequently exposed HiFi
challenges and 24 dependent depth/name perturbations are development evidence.
The six challenges are no longer held out. The following table re-scores archived
before/after outputs with one strict evaluator; it does not imply every current
source revision was rerun over every sample.

| Cohort | Inputs | Exact individual sequences, before → after | Exact diploid sequences | Exact mutation annotations, before → after | Extra mutation records, before → after |
| --- | ---: | --- | --- | --- | --- |
| Original HiFi | 44 | 39/88 → 39/88 | 10/44 → 10/44 | 24/28 → 24/28 | 14 → 14 |
| Original ONT | 3 | 3/6 → 3/6 | 0/3 → 0/3 | 2/2 → 2/2 | 0 → 0 |
| Later HiFi challenges | 6 | 1/12 → 1/12 | 0/6 → 0/6 | 2/3 → 2/3 | 0 → 0 |
| Dependent perturbations | 24 | 6/48 → 7/48 | 0/24 → 0/24 | 6/14 → 7/14 | 12 → 11 |

The perturbation improvement includes repairing an invalid legacy alias artifact;
it is not proof of a new reconstruction algorithm. The six low-depth no-calls
remain in denominators. Matching is one-to-one and retains assignment ambiguity;
an unproven duplicate sequence receives no independent second-allele credit.

Historically, 15/16 HiFi normal samples had no mutation alarm. Ten have unresolved
reconstruction under the stricter contract. This leaves five confident negatives,
one false-alarm control and ten unresolved negatives. Conditional specificity is
5/6, with the unresolved fraction reported separately. The single ONT normal
control is unresolved, so its conditional specificity is undefined. These small,
simulation-based cohorts do not establish clinical sensitivity or specificity.

## Experiments that were not promoted

| Experiment | Observed result | Decision |
| --- | --- | --- |
| Minimum coverage 10 → 5 | More counts and outputs, but extra events 25 → 38 and normal false-alarm samples 4/30 → 7/30 across 77 inputs | Retain default 10 |
| Optional read-backed phasing | Cached exact individual sequences 1/10 → 4/10, but an extra false event; no exact diploid recovery | Experimental opt-in, disabled by default |
| Minimum selected-phase support | Can remove a false singleton but also reject real minority alleles; original 30:2 became internally selected 14:1 | No universal support floor |
| Strict segmentation stop | Lost five named events among 106 exposed consensus outputs | Experimental only; retain default recovery and expose unresolved sequence |
| Simple anchor length modes | Promising wider-gap counts, but bin suppression merges true one-repeat-separated alleles; many ONT reads lack usable anchors | Diagnostic only |
| Continuous-span one/two component models | Improved count-only results but still split the equal-length control into 59/60; downstream mutation recovery untested | Not promoted |

The equal-length failure has strong reference-fit artifacts: the false longer
candidate has many alignment records but little full-boundary read support. The
25/140 case loses all five long primary records before calling because their
reference-fit counts are spread below the per-contig threshold; its false second
short component is fed one primary read. These explain the failures without
establishing a safe replacement inference/assignment method. The tests retain
both strict expected failures and their original tolerances.

## Reproducing and extending the evidence

Use the [simulation experiment guide](simulation-experiments.md) to specify new
HiFi and ONT cases, seeds, lengths, mutation targets and requested templates in
JSON. Keep simulator and platform model configurations explicit. Requested
coverage is not the number of retained reads: record actual usable inputs.

Freeze caller settings, tool/model versions, source revision and evaluation rules
before creating a final panel. Keep generation truth out of calling logic. If you
change settings after examining final results, classify that panel as development
and use new seeds for the next final evaluation. Report sample failures, no-calls,
assignment ambiguity, extra calls and per-sample regressions alongside aggregates.
Perturbations from the same source sample are dependent observations.

The repository's `.planning/` evidence reports retain experiment decisions,
commands, input/source hashes and reviewer dispositions. Raw reads, results,
indexes and model files remain local generated artifacts. See also
[benchmarking](benchmarking.md), [configuration](configuration.md), and
[limitations](../reference/limitations.md).
