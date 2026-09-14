# Read-backed phasing feasibility evidence — 2026-09-14

Bounded investigation requested before final source freeze. Production source and
tracked tests were not edited. All executable prototypes, fixtures and reports
are ignored under `tests/results/production_validation_20260914/read_phasing/`.
The earlier premise that phasing tooling was unavailable is false.

## Actual installed capability and commands

`/home/bernt-popp/miniforge3/envs/env_clair3/bin/whatshap --version` returned
**1.7**. `whatshap phase --help` was captured in `help.txt`; default mapping quality
is 20, default maximum internal coverage is 15, and indels are excluded unless
`--indels` is supplied. The environment also supplies Clair3's final-phasing flag,
but simply enabling that flag would not address colliding record names.

Probes invoked the existing `muc_one_span.tools.run_tool` boundary with argument
lists, including:

```
whatshap phase --indels --ignore-read-groups --sample sample \
  --reference REFERENCE --output-read-list READ_LIST -o PHASED_VCF INPUT_VCF BAM
```

No genotype distrust, homozygous inclusion, read merging, lowered MAPQ or altered
upstream remapping was enabled. Synthetic exact indel sequences were generated
using actual `bcftools consensus -s sample -H 1/2`; names were compared as an
unordered diploid pair because swapping both haplotype labels is immaterial.

Executed commands, all successful after explicitly recording the empty-BAM tool
failure in the synthetic report:

- `uv run --locked --all-extras python tests/results/production_validation_20260914/read_phasing/probe.py`
- `uv run --locked --all-extras python tests/results/production_validation_20260914/read_phasing/indel_probe.py`
- `uv run --locked --all-extras python tests/results/production_validation_20260914/read_phasing/cached_probe.py`

## Record identity contract and sensitivity

The prototype reads existing coordinate-sorted BAM records through samtools,
keeps mapped primary records, replaces each QNAME by `phase_record_<ordinal>`,
and writes a JSON source map. Every non-QNAME SAM field is asserted identical
before and after conversion; this includes sequence, quality, CIGAR, coordinates,
flags and tags. Indexing is the only downstream alignment operation. The current
prototype source ordinal counts the SAM stream including headers; production
should explicitly count all alignment records and name that field accordingly.
The source-map record order still identifies every retained primary record.

The duplicate-QNAME synthetic fixture deliberately uses the identical name for
20 independent primary records. Without renaming, WhatsHap rejects it because
that name occurs more than twice. Renaming preserves all 20 records and produces
correct cis/trans phase. Source-haplotype lookup on the tool's read list gives
pure assignments to both haplotypes for cis, trans and linked partial reads.
The list contains 15 selected reads in the complete two-site fixtures due to
WhatsHap's default internal downsampling; 20 retained BAM records must not be
reported as 20 independently used phasing observations.

A paired unique-name versus renamed-name experiment for cis/trans/shared indels
recovers exactly the same unordered two full sequences in all three cases.
Whole-pair haplotype labels can swap. Thus literal VCF GT-string equality would
be an incorrect invariance test. This bounded result does not establish arbitrary
high-depth/noisy VNTR name-order insensitivity.

## Synthetic results

| Fixture | Output phase contract | Sequence / assignment result |
| --- | --- | --- |
| Two SNPs cis | Common PS, phased | Exact unordered sequence pair; source-pure assignments |
| Two SNPs trans | Common PS, phased | Exact unordered sequence pair; source-pure assignments |
| Shared homozygous SNP plus heterozygous SNP | Single heterozygous unordered | Exact pair; no multi-site phase asserted |
| Single multiallelic 1/2 | Single heterozygous unordered | Exact pair; no reference allele substituted |
| Multiallelic 1/2 plus second heterozygous site | Unphased | **Unsupported by installed tool**, no haplotype assertion |
| Unphased VCF, zero BAM records | Explicit tool failure | Production should preflight and retain unresolved evidence |
| Reads cover sites separately, no linking read | Unphased | No multi-site phase assertion |
| Partial reads spanning both sites | Common PS, phased | Exact pair; source-pure assignments |
| Partial GT 0/. | Missing genotype | No haplotype assertion |
| Empty VCF | No informative heterozygosity | No sequence identity or homozygosity assertion |
| Two insertions cis | Common PS, phased | Both bcftools sequences exact, original and renamed names |
| Two insertions trans | Common PS, phased | Both bcftools sequences exact, original and renamed names |
| Shared insertion plus heterozygous insertion | Single heterozygous unordered | Both bcftools sequences exact, original and renamed names |

`synthetic_report.json` includes a mechanical GT-index replay equality field for
all fixtures. Equality in an unresolved fixture is merely incidental input GT
ordering and **must not be interpreted as recovered phase or an exact caller
output**. Source assignment purity is meaningful only for the multi-site phased
fixtures; one-site and empty cases do not have a two-group phasing read list.
These fixtures use error-free random-sequence alignments, not repetitive noisy
VNTRs. They prove software feasibility, not general biological performance.

## Cached development candidates

Full input paths, original/renamed VCF records, commands and read counts are in
`cached_report.json`. No cached truth read-source map exists, so assignment purity
cannot be inferred from QNAME, sequence similarity or the output haplotype tag.

| Cached candidate | Primary records / original unique names | Before → after renamed phasing |
| --- | --- | --- |
| normal_60_80 n20 seed1703 | 12 / 12 | Unphased → phased |
| dupc_60_80 n20 seed1701 | 13 / 13 | Unphased → phased |
| normal_60_80 n20 seed1701 | 11 / 11 | Unphased → phased |
| dupc_60_80 n20 seed1702 | 11 / 11 | No heterozygosity → unchanged |
| heldout_normal_60_62 (already exposed development) | 180 / 102 | Unphased → phased |
| canonical homozygous60_60 allele_1 | 164 / 99 | Unphased → phased (diagnostic only) |
| canonical homozygous60_60 allele_2 | 2 / 2 | Unphased → unchanged (diagnostic only) |

Original-name and renamed-name cached runs all completed, but were not scored
against exact haplotype truth in this bounded probe. Phase acquisition alone is
not an accuracy result. Existing false single-length assignment, genotype errors,
repeat-alignment artifacts and reference-fill errors remain possible.

Canonical60_60 currently has `same_length=false`, inferred lengths 59 and 62 and
contigs 51/54. The proposed D hook would never execute for that input. Diagnostic
phase of its separate BAMs cannot repair the upstream false split or justify
emitting more alleles. This limitation must remain visible.

## Proposed production boundary and promotion decision

A small dedicated read-phasing module can be called **after filtered variants are
available and only within same-length disambiguation** when multiple unresolved
heterozygous records are present. Preserve already supported common-PS records,
no-variant, single-heterozygous, missing/non-diploid/conflicting records. Select the
VCF sample explicitly, create collision-free primary BAM names plus source map,
preflight zero usable records, invoke the tool using the existing abstraction,
compress/index its output using existing bcftools, and then apply the unchanged
`phase_evidence` guard. Persist tool version, original/selected primary record
counts, source-map path, phased VCF path and read-list path in additive evidence.
Never interpret a PS as reference confidence or molecular independence.

The tool is installed transitively here, but production availability must be
explicitly checked and documented. If not available, retain the original VCF and
an explicit unavailable status; actual subprocess failures must remain visible.
Do not silently introduce a new dependency or change Clair3 mapping/calling.
Empty BAM is insufficient evidence and should be handled before invoking a tool
known to reject it. A wrapper should avoid overwriting existing VCFs and BAMs.

**Feasible for biallelic linked SNPs/indels, with explicit unsupported multiallelic
multi-site status. Not yet promoted.** To promote: root must authorize source
ownership, add failing deterministic wrapper/record tests, real integration tests,
and exact cached/fresh evaluation of the hooked production path. The present
probe alone cannot meet the broad no-regression scientific promotion gate; it
removes the tooling objection and supplies concrete implementation evidence.
