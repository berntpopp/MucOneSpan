# MucOneSpan

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/berntpopp/MucOneSpan/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/MucOneSpan/actions)
[![Documentation](https://img.shields.io/badge/docs-MkDocs%20Material-blue)](https://berntpopp.github.io/MucOneSpan/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Open-source MUC1 VNTR analysis pipeline for PacBio HiFi and ONT amplicon data.

---

## Overview

MucOneSpan analyzes long-read sequencing data of the MUC1 VNTR region. It implements and extends the methods described by Vrbacka et al. (2025), with support for PacBio HiFi and ONT amplicons. [MucOneUp](https://github.com/berntpopp/MucOneUp) supplies simulated haplotypes and reads for benchmarking.

**Key Capabilities:**

- Detect frameshift mutations in the MUC1 VNTR that cause ADTKD-MUC1 kidney disease
- Resolve close allele pairs (3-9 repeats apart) via indel-valley splitting
- Pre-computed mutation template catalog (13 known mutations) with O(1) lookup
- Per-repeat and per-allele confidence scoring with VCF cross-validation
- Fully automated -- no manual IGV inspection steps

---

## Quick Start

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi \
  --threads 8
```

**Documentation:** https://berntpopp.github.io/MucOneSpan/

---

## Installation

Install the Python package with HTML report support:

```bash
pip install 'muc_one_span[report] @ git+https://github.com/berntpopp/MucOneSpan.git@v0.10.0'
muconespan --version
```

The display name is **MucOneSpan**, the command is `muconespan`, and the Python
package/import is `muc_one_span` (normalized distribution name: `muc-one-span`).
The current release is **0.10.0**.

For development:

```bash
git clone https://github.com/berntpopp/MucOneSpan.git
cd MucOneSpan
make dev
uv run muconespan --version
```

**Requires:** minimap2, samtools, bcftools, Clair3 on PATH. See [Installation Guide](https://berntpopp.github.io/MucOneSpan/getting-started/installation/).

---

## Citation

If you use MucOneSpan, please cite the software and the original method:

**Software:**

```bibtex
@software{popp_muc_one_span_2026,
  author       = {Popp, Bernt},
  title        = {MucOneSpan: Open-source MUC1 VNTR analysis pipeline},
  version      = {0.10.0},
  year         = {2026},
  url          = {https://github.com/berntpopp/MucOneSpan},
}
```

**Original method:**

> Vrbacka A, Pristoupilova A, Kidd KO, et al. Long-Read Sequencing of the MUC1 VNTR: Genomic Variation, Mutational Landscape, and Its Impact on ADTKD Diagnosis and Progression. *bioRxiv.* 2025. doi: [10.1101/2025.09.06.673538](https://doi.org/10.1101/2025.09.06.673538)

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

---

## License

MIT License -- see [LICENSE](LICENSE) for details.

---

**Maintained by:** [Bernt Popp](https://github.com/berntpopp)
