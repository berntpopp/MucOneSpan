# Milestone 1 consensus-command equivalence

Date: 2026-09-14

## Question and method

The milestone review asserted that changing the ordinary per-length consensus
command from:

```text
bcftools consensus -f REFERENCE VCF
```

to:

```text
bcftools consensus -s SAMPLE -H I -M N -f REFERENCE VCF
```

would stop applying heterozygous indels. This was tested directly with installed
bcftools 1.17, holding each cached filtered `variants.vcf.gz` and matching
single-contig `ref_contig_*.fa` fixed. The matrix was the complete original
development set under `tests/results/deep_validation_20260914/{hifi,ont}`: 44
HiFi and three ONT samples, two allele VCF/reference pairs per sample. No mapping,
variant calling, filtering, trimming, classification, or caller algorithm ran.

The ignored driver and machine-readable output are:

- `tests/results/deep_validation_20260914/consensus_equivalence_20260914/run.py`
- `tests/results/deep_validation_20260914/consensus_equivalence_20260914/results.json`
- results JSON SHA-256:
  `a4cfec59aa536cb4d8e8a95b0012cc276c2a1ef2a4c00d4006cb1570c1474964`

Each result row records input reference/VCF hashes, selected VCF sample, both
output hashes, byte and full-sequence equality, sequence lengths, first changed
position, ambiguous-base counts, and heterozygous-indel records.

## Results

| Platform | Samples | Alleles | Byte-identical | Full-sequence identical | Het-indel records | Alleles with het indels |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HiFi | 44 | 88 | 88 | 88 | 17 | 9 |
| ONT | 3 | 6 | 6 | 6 | 0 | 0 |
| Total | 47 | 94 | 94 | 94 | 17 | 9 |

There are no changed cases. `changed_cases` is an empty array. The full FASTA
byte comparison is stronger than sequence equality and therefore also rules out
header or line-wrapping changes for this matrix.

All 17 heterozygous indels are genotype `0/1`. They occur in:

| Sample | Allele | Records | Positions and alleles |
| --- | --- | ---: | --- |
| `sample_bench_5003` | allele_2 | 1 | 7392 GC>G |
| `sample_bench_5004` | allele_1 | 1 | 1392 G>GC |
| `sample_dupc_40_50` | allele_1 | 1 | 1392 G>GC |
| `sample_dupc_50_55` | allele_2 | 3 | 1092, 1692, 2412 GC>G |
| `sample_dupc_60_80_s3` | allele_2 | 1 | 4692 GC>G |
| `sample_dupcccc_60_80` | allele_2 | 1 | 4692 GC>G |
| `sample_gap4_40_44` | allele_2 | 5 | 1032, 2412, 2712, 2892 GC>G; 2964 TCCAC>T |
| `sample_insg_100_120` | allele_2 | 3 | 4092, 5652, 5712 GC>G |
| `sample_short_25_30` | allele_1 | 1 | 1797 GC>G |

## Interpretation

The review's specific prediction does not reproduce: bcftools 1.17 `-H I`
applies these `0/1` indels exactly as the historical no-sample/no-haplotype
command does. This agrees with the separately executed minimal `0/1` indel
fixture, in which both commands selected ALT.

The result establishes command equivalence for the cached single-sample,
QUAL-filtered 44+3 matrix. It does not generalize to multi-sample VCFs, missing
genotypes, other bcftools versions, or genotype patterns absent from these 94
VCFs. Explicit sample selection still improves determinism and rejects ambiguous
multi-sample inputs. No calling or consensus algorithm change follows from this
experiment.
