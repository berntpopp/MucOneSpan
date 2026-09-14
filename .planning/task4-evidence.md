# Task 4 — consensus and phase evidence

Implemented in the production-validation worktree, 2026-09-14. No commits,
dependency changes, threshold changes, or read-name transformation.

## Behavior and interfaces

- VCF selected-sample parsing now retains CHROM, POS, REF, comma-separated ALT,
  QUAL (`None` when missing), GT, PS (`None` when absent), and sample. The GT
  separator remains in `genotype`. Query failures and malformed records raise;
  multi-sample input needs explicit selection. QUAL/PASS filtering is unchanged.
- `phase_evidence` requires complete diploid genotypes and one chromosome/PS
  block across informative heterozygous records before selecting two ordered
  genotype haplotypes. One heterozygous site yields an unordered pair. Missing,
  haploid, unphased, disconnected, missing PS and overlapping records have
  explicit states. Overlapping reference intervals are conservatively unresolved.
- Equal-length calling no longer creates WT/all-ALT VCFs. A supported pair
  shares its filtered VCF and selects actual GT positions 1/2. Other evidence
  yields one mixed candidate. Empty and all-shared calls do not prove sequence
  homozygosity. A legacy second allele remains marked `not_separately_resolved`
  and is absent from the returned VCF map. The path map remains path-only.
- Ordinary per-length calls use genotype-aware IUPAC and retain residual
  heterozygosity status. `build_consensus(..., sample=None, haplotype='I')`
  explicitly selects sample/GT policy; missing genotypes are masked with N.
  IUPAC cannot express heterozygous indels; these remain candidate artifacts,
  not claims of resolved indel phase.
- All produced candidates have `sequence_identity_status='unresolved'`,
  `reference_confidence='unverified'`, and a candidate reconstruction status.
  `variant_observation_group` is the source VCF; `sequence_source` appends
  `:GT1`, `:GT2`, or `:GTI`. `independent_haplotype_evidence` is true only for
  genotype-supported selected pair members. It does not claim molecule/read-set
  independence. Two sources from a shared VCF remain one variant observation group.
- `allele_info['consensus_context']` contains `vcf_path`, `reference_path`,
  `full_consensus_path`, `chrom`, `sample`, `haplotype`, and actual full-consensus
  `trim_start`/`trim_end` (0-based, half-open). The same dict is written to
  `consensus_<allele>_context.json`. Paths identify runtime artifacts. This
  context enables the coordinator's replay-verified projection engine; the
  sidecar alone does not establish a valid projection.

## Observed red tests

Before implementation, focused changed VCF/calling/consensus tests had 14 failures
and 44 passes: swallowed query errors, malformed-data skipping, discarded identity,
implicit sample/haplotype, empty-call homozygosity and unsafe split metadata.
Nine same-length phase scenarios then failed against the old implementation.
Two tests failed for implicit IUPAC and absent actual trim intervals; one failed
for absent ordinary residual-heterozygosity state; one exposed overlapping
records incorrectly labeled phased. New provenance assertions failed before
sidecar/source metadata was added. Temporary logs: `/tmp/task4-*-red.log`.

## Verification

Executed with real bcftools 1.17 and samtools on the env_clair3 PATH, preserving
its Python for external tools:

```bash
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH \
uv run --locked --all-extras pytest \
  tests/unit/test_vcf.py tests/unit/test_calling.py \
  tests/unit/test_consensus.py tests/unit/test_phasing.py \
  tests/integration/test_phase_consensus.py tests/integration/test_consensus.py \
  --no-cov -q
```

Result: **86 passed**, including **71 unit and 15 real-tool integration tests**,
zero skips. Fixtures cover exact trans/cis, shared variants, GT1|2, 2/2, 0/0,
single unphased heterozygous insertion, missing GT, haploid GT, no records,
unphased/missing/disconnected PS, multi-sample rejection/selection, and context
round-trip. Unit tests cover malformed/query failures and overlapping records.

Owned files pass `uv run --locked --all-extras ruff check ...`,
`ruff format --check ...`, and configured `mypy` on the four changed/new package
modules. All owned code files are under 650 lines (largest: test_calling.py, 480).
A concurrent broad-unit snapshot had 315 passes/9 failures across changing
classification, evaluation and CLI tests; the coordinator was notified. That
snapshot is not a final global check. Coordinator owns make ci-check, full
integration/build gates and production-pipeline evaluation.

## Review finding

The plan review's claim that bcftools `-H1` errors for unphased `0/1` is rejected
by an actual installed-tool fixture: reference AAAAAAAAAA, POS2 A>AC, GT0/1
produces AAAAAAAAAA with `-s SAMPLE -H1` and AACAAAAAAAA with `-H2`. Those two
outputs are an unordered pair at one heterozygous site. Multiple unphased sites
still cannot establish cross-site phase.

## Safety changes and unmet science

The old fabricated two-allele WT/all-ALT outputs are replaced by one uncertain
candidate when multiple loci are unphased, missing, disconnected or conflicting.
Missing and haploid genotypes also do not establish a diploid pair. Empty/all-hom
calls retain one candidate with unknown sequence identity. These reduce false
assertions and can increase missing-allele/no-call counts; they are safety-contract
changes, not demonstrated scientific accuracy gains. The evaluator must keep
these samples and all expected alleles in its denominators.

Read-backed phase is **not implemented or promoted**. This change only interprets
existing selected-sample VCF evidence. Existing extraction/remapping does not
provide validated collision-free molecule identity or a tested read-assignment
phase graph; QNAME deduplication would lose distinct source records in known
simulations. Neither software fixtures nor retained PS tags establish read-backed
phase sensitivity on HiFi/ONT VNTRs. Whole-reference confidence, uncovered-base
masking outside explicit missing GT, and clinical calibration remain unverified.
