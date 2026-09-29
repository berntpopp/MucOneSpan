# MucOneSpan

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/berntpopp/MucOneSpan/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/MucOneSpan/actions)
[![Documentation](https://img.shields.io/badge/docs-MkDocs%20Material-blue)](https://berntpopp.github.io/MucOneSpan/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MucOneSpan analyzes the MUC1 VNTR from PacBio HiFi and Oxford Nanopore amplicon reads, reporting allele lengths, repeat composition, frameshift variants, and confidence scores.

## Install

```bash
pip install 'muc_one_span[report] @ git+https://github.com/berntpopp/MucOneSpan.git@v0.17.0'
```

The default hybrid engine needs no external tool for FASTQ input (`samtools` for BAM input). Its POA library `pyabpoa` builds from source and needs a C compiler and zlib. See the [installation guide](https://berntpopp.github.io/MucOneSpan/getting-started/installation/) for setup and container options.

## Run

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --report
```

Since 0.17.0, `muconespan run` uses the read-centric **hybrid engine** by default for PacBio HiFi and ONT reads, amplicon and genomic input (`--assay` is recorded for provenance only). The hybrid engine needs no platform setting and runs single-threaded; the ladder-only options (`--platform`, `--threads`, `--min-coverage`, `--mapping-timeout`, `--reference`, `--clair3-model`, `--min-qual`, `--minimap2-preset`) are ignored with a warning and recorded in `summary.json["ignored_options"]`. For a WGS BAM, pass a BAM already subset to the MUC1 region: the hybrid engine streams every read of the BAM.

The ladder engine (minimap2, Clair3 and bcftools) is **deprecated**: `--engine ladder` still works, prints a warning and records it in `summary.json`; it will be removed in a later release. See the [migration guide](https://berntpopp.github.io/MucOneSpan/guides/migration/) and the [configuration guide](https://berntpopp.github.io/MucOneSpan/guides/configuration/#hybrid-engine).

## Acknowledgment

MucOneSpan draws on some ideas from [Vrbacka et al. (2025)](https://doi.org/10.1101/2025.09.06.673538).

[Documentation](https://berntpopp.github.io/MucOneSpan/) · [Development](docs/development.md) · [MIT license](LICENSE)


Design reproducible HiFi and ONT simulation panels with the
[experiment guide](docs/guides/simulation-experiments.md). Read the
[0.11.0 validation results](docs/guides/validation-results.md) and the
[clinical validation results](docs/guides/clinical-validation-results.md) for measured
real-world performance on PRJEB92208 data and comparator evaluation.
