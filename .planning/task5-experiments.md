# Task 5 E/F: diagnostic spans and local indel evidence

Date: 2026-09-14. Exposed development data only. No length, QUAL, calling,
consensus, or assignment defaults changed by this experiment. No rescue promoted.

## Diagnostic implementation and contract

`src/muc_one_span/read_evidence.py` provides pure `spanning_evidence` and
`assign_span` functions. It does not read simulator truth, files, BAMs, or execute
tools. Callers supply original source records exactly once; when supplying BAM
records, exclude secondary and supplementary copies. Record index, not QNAME,
identifies an observed input record. Distinct records with identical names and
sequences remain distinct. This does not establish independent PCR molecules.

Anchors use literal Levenshtein substring distance with a specified maximum
per-anchor edit allowance and both read orientations. Only minimum-edit matches
for each anchor are retained. Multiple distant loci or orientations are ambiguous;
nearby tied optimal boundaries produce an explicit span interval. Missing one
anchor, missing both, and ambiguity are explicit rejection reasons. The partial
status includes anchors that cannot form an ordered nonoverlapping pair; it is
not proof of physical truncation. A degraded true anchor can lose to a better
ectopic match. Neither error tolerance nor interval width is a calibrated
probability, and interior sequencing indels are not corrected.

Span assignment returns all candidate lengths compatible with the interval and
an explicit base-distance allowance. Equal lengths remain ambiguous between
haplotypes. A singleton candidate does not establish sequence homozygosity.
The functions are diagnostic APIs and are not wired into production inference.

## Retained input matrix and matching experiment

Executed all 77 retained input runs: original 44 HiFi, original 3 ONT, previously
exposed 6 HiFi, and 24 prior perturbations. The latter reuse reads and are not
independent simulations. Source BAMs were streamed with `samtools view -F 2304`,
retaining primary unmapped records; FASTQs were read by record. No QNAME collapse.
Input SHA-256 values and individual record evidence are in ignored artifacts.
Original source assignment truth is unavailable; assignment purity, wrong-source
rates, and haplotype-specific support cannot be measured on these 77 runs.

Anchors are the first 20 bases of dictionary repeat 1 and last 20 bases of repeat
9. Compared exact matching with allowances of one and two edits per anchor.
Inference runs before truth loading. Only afterward are distinct-length candidate
sets compared to stored simulation counts. Rounded span/60 is exploratory;
mutated unit length and sequencing indels mean it is not sequence reconstruction.

| Cohort | Records | Exact spanning | ≤1 edit spanning | ≤2 edits spanning |
| --- | ---: | ---: | ---: | ---: |
| Original HiFi, 44 | 7,257 | 6,241 (86.0%) | 6,906 (95.2%) | 6,986 (96.3%) |
| Original ONT, 3 | 600 | 112 (18.7%) | 209 (34.8%) | 321 (53.5%) |
| Exposed HiFi, 6 | 1,033 | 885 (85.7%) | 992 (96.0%) | 1,000 (96.8%) |
| Perturbations, 24 | 1,403 | 1,213 (86.5%) | 1,331 (94.9%) | 1,346 (95.9%) |

No ambiguous anchor-pair reads arose in these inputs; deterministic fixtures
exercise both repeated loci and conflicting orientations. At two edits, 157,
63, 24, and 31 spanning records respectively have tied boundary uncertainty.
ONT still rejects 279/600 (46.5%), so tolerance does not solve coverage bias.

Candidate generation chooses highest-support exact length bins, ties favoring
shorter lengths, and suppresses bins within ±gap. It keeps all supported modes,
without forcing two alleles. All combinations of support 3/5/10, edit allowance
0/1/2, and suppression gap 0/1/2/3 were retained, without selecting new defaults.
The following fixed support-five comparison illustrates the tradeoff:

| Cohort/matching | gap 0 exact sets | gap 1 | gap 2 | gap 3 |
| --- | ---: | ---: | ---: | ---: |
| Original HiFi, exact | 30/44 | 43/44 | 43/44 | 41/44 |
| Original HiFi, ≤2 edits | 23/44 | 43/44 | 43/44 | 41/44 |
| Original ONT, exact | 2/3 | 3/3 | 3/3 | 3/3 |
| Original ONT, ≤2 edits | 0/3 | 3/3 | 3/3 | 3/3 |
| Exposed HiFi, exact | 3/6 | 5/6 | 4/6 | 3/6 |
| Exposed HiFi, ≤2 edits | 2/6 | 5/6 | 4/6 | 3/6 |
| Perturbations, exact | 12/24 | 13/24 | 13/24 | 13/24 |
| Perturbations, ≤2 edits | 12/24 | 14/24 | 14/24 | 14/24 |

Gap zero admits error bins: at two edits it adds 27 extra HiFi lengths, six ONT
lengths, four exposed lengths, and two perturbation lengths. Gap one misses the
exposed 60/61 distinction; gap two also misses 60/62; gap three also misses
60/63. More anchoring does not establish the correct number of alleles. The
support-five original HiFi miss remains the weak 140-repeat component of 25/140.

At support five/gap one with ±30-base compatibility, unique compatible records
rise from 5,995 to 6,636 HiFi, 89 to 219 ONT, 774 to 866 exposed, and 1,145 to
1,261 perturbation records. There are no multi-candidate compatible records
under that separated-mode rule. These are candidate compatibility counts, **not
assignment accuracy**; the missing candidate can leave its true reads unassigned
or compatible with another length. Scoring uses distinct length sets, not a
claim to distinguish the two sources of a same-length pair.

## Known-source deterministic fixtures

Small exact synthetic strings and the unit suite cover true gaps 0/1/2/3,
colliding names, equal-length distinct sequences, interior deletion/insertion,
anchor substitution/insertion/deletion, both orientations, partial reads,
multiple loci, and assignment ties. The independent substring DP oracle checks
the bit-vector anchor search exhaustively over 64 short strings, three anchors,
and edit allowances one and two.

The ignored known-source experiment uses 12:12 and 30:3 source-record mixtures
with full, anchor-insertion, interior-deletion, reverse and two terminal-partial
read types. Every record has the same name. At gap zero, 16/24 or 23/33 records
are ambiguous and the rest unassigned. At gaps 1/2/3, those 16 or 23 are correctly
compatible with their known source, zero wrong assignments, and the same partial
reads remain unassigned. This tests contract behavior using explicit candidate
spans; it does not establish source recovery or calibration on simulated HiFi/ONT.

## Candidate-local CIGAR investigation

Read actual existing remapped allele BAMs for two positives and HiFi/ONT normal
controls, retaining each primary mapped alignment record. The experiment inspects
only predeclared loci and their ±90-base neighborhoods. It left-normalizes CIGAR
insertion sequences against the actual single-contig reference. Its coordinates
are interbase offsets, numerically equal to the preceding VCF anchor POS for
these insertions. No motif-anywhere scan or universal QUAL change is used.

| Input/allele | Candidate VCF anchor | Single-C insertion records | Total primary records |
| --- | ---: | ---: | ---: |
| dupC100/120 allele 1 | contig_91:3192 | 84 | 111 |
| dupC100/120 allele 2 | contig_111:3192 | 2 | 56 |
| long120/140 allele 1 | contig_111:3492 | 43 | 110 |
| long120/140 allele 2 | contig_131:3492 | 1 | 64 |
| HiFi normal60/80 allele 1 | contig_51:3192 / 3492 | 0 / 0 | 108 |
| HiFi normal60/80 allele 2 | contig_71:3192 / 3492 | 3 / 4 | 55 |
| ONT normal60/80 allele 1 | contig_51:3192 / 3492 | 0 / 0 | 127 |
| ONT normal60/80 allele 2 | contig_69:3192 / 3492 | 2 / 1 | 61 |
| HiFi normal60/60 allele 1 | contig_51:3192 / 3492 | 10 / 10 | 164 |

The existing 100/120 raw record is `G>GC`, QUAL 3.22, GT1/1, DP110, AD22,84.
All 84 insertion-supporting alignments have MAPQ60. One of the 111 mapped records
has MAPQ1, explaining why the raw-record denominator need not equal caller DP.
The raw candidate is real and filtered by the existing QUAL5 threshold.

The 120/140 missing locus has 43 single-C plus five double-C insertions at 3492,
25 single-C insertions at adjacent repeat anchor3432, and eight at3552. At3492,
60 records have a continuous local match block, six other spanning configurations,
and one deletion over the locus. A match block is aligned coverage, **not proof
of the REF sequence**. The pileup VCF contains a `RefCall` at3492 with QUAL15.65,
GT0/0, DP110, AD61 and reference fraction0.5545. There is no full-alignment record
at this locus and no final indel. Therefore the failure occurs despite read-level
insertion evidence and includes candidate/model selection, rather than a blanket
absence of evidence in long reads. Most alignments have MAPQ60, but single-contig
MAPQ does not establish unambiguous position within similar repeat units.

A targeted variant rescue remains **unmet**. The mapped positive locus is informed
by the prior truth audit, not a validated caller-side candidate-discovery rule.
Normal controls have nonzero equivalent insertion support, nearby repeat anchors
compete, original source assignment is missing, and a CIGAR observation alone
cannot establish exact repeat localization or genotype. No rescue VCF, consensus,
mutation label, sensitivity gain, or new cutoff is asserted. Further work needs
candidate REF/ALT local sequence likelihoods, repeat-placement ambiguity and
source-aware controls; neither lowering QUAL globally nor interpreting all CIGAR
insertions as known mutations meets that requirement.

## Artifacts and executed validation

Ignored root: `tests/results/production_validation_20260914/read_evidence/`.
`anchor_experiment.py`, `anchors.json`, `anchors_summary.json`, per-input evidence,
`known_source_experiment.py`, `known_source.json`, `cigar_experiment.py`,
`cigar_summary.json`, record-level CIGAR projections and logs retain the work.
`source_hashes.json` fingerprints source and experiment scripts. Absolute local
paths occur only in ignored scripts/artifacts, not authored package or tests.

Commands actually executed:

- `uv run --locked --all-extras pytest tests/unit/test_read_evidence.py --no-cov -q`:
  initial missing-module collection failure, then 12, 14, and 19 passing tests
  as implementation and additional cases were added.
- `uv run --locked --all-extras pytest -o addopts='' tests/unit/test_read_evidence.py --cov=muc_one_span.read_evidence --cov-branch --cov-report=term-missing --cov-fail-under=80 -q`:
  final **20 passed**, **100%** statement/branch coverage for this module.
- A prior focused coverage command retained repository-wide coverage addopts:
  all 19 tests passed but aggregate coverage failed at9%; this was a command-scope
  mistake, not a passing full-suite gate. No aggregate gate was weakened.
- Focused `ruff check`, `ruff format --check`, and `mypy` passed for owned files.
- Three ignored experiment scripts executed with locked all-extras Python. The
  anchor sweep covered all77 inputs in about89 seconds of measured sample work.
  CIGAR inspection exercised installed samtools1.15.1 through `tools.py`.

Coordinator owns the final `make ci-check`, integration/build checks, and final
source-wide diff review. This worker has not committed or changed dependencies.
