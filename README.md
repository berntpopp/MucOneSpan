# MucOneSpan

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/berntpopp/MucOneSpan/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/MucOneSpan/actions)
[![Documentation](https://img.shields.io/badge/docs-MkDocs%20Material-blue)](https://berntpopp.github.io/MucOneSpan/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MucOneSpan analyzes the MUC1 VNTR from PacBio HiFi and Oxford Nanopore amplicon reads, reporting allele lengths, repeat composition, frameshift variants, and confidence scores.

## Install

```bash
pip install 'muc_one_span[report] @ git+https://github.com/berntpopp/MucOneSpan.git@v0.13.0'
```

The pipeline also requires minimap2, samtools, bcftools, and Clair3. See the [installation guide](https://berntpopp.github.io/MucOneSpan/getting-started/installation/) for setup and container options.

## Run

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi \
  --threads 8 \
  --report
```

For Oxford Nanopore reads, add `--platform ont` and use a matching Clair3 model.

## Acknowledgment

MucOneSpan draws on some ideas from [Vrbacka et al. (2025)](https://doi.org/10.1101/2025.09.06.673538).

[Documentation](https://berntpopp.github.io/MucOneSpan/) · [Development](docs/development.md) · [MIT license](LICENSE)


Design reproducible HiFi and ONT simulation panels with the
[experiment guide](docs/guides/simulation-experiments.md). Read the
[0.11.0 validation results](docs/guides/validation-results.md) for measured
improvements, rejected experiments and remaining reconstruction limitations.
