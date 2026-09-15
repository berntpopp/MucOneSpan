# Citation Guide

How to cite MucOneSpan in your research publications.

---

## Acknowledgment

MucOneSpan draws on some ideas from Vrbacka et al. (2025):

### Vrbacka et al. 2025

```bibtex
@article{vrbacka2025muc1,
  author = {Vrbacka, A. and Pristoupilova, A. and Kidd, K.O. and others},
  title = {Long-Read Sequencing of the MUC1 VNTR: Genomic Variation,
           Mutational Landscape, and Its Impact on ADTKD Diagnosis
           and Progression},
  journal = {bioRxiv},
  year = {2025},
  doi = {10.1101/2025.09.06.673538}
}
```

> Vrbacka A, Pristoupilova A, Kidd KO, et al. Long-Read Sequencing of the MUC1 VNTR: Genomic Variation, Mutational Landscape, and Its Impact on ADTKD Diagnosis and Progression. *bioRxiv.* 2025. doi: [10.1101/2025.09.06.673538](https://doi.org/10.1101/2025.09.06.673538)

---

## Software Citation

Software citation information:

```bibtex
@software{muconespan2026,
  author = {Popp, Bernt},
  title = {MucOneSpan: Open-Source MUC1 VNTR Analysis Pipeline},
  version = {0.13.0},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/berntpopp/MucOneSpan},
  note = {Software version available at https://github.com/berntpopp/MucOneSpan/releases}
}
```

---

## MucOneUp Citation

If you use MucOneUp for benchmarking or test data generation:

```bibtex
@software{muconeup2025,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/berntpopp/MucOneUp},
  note = {Software version available at https://github.com/berntpopp/MucOneUp/releases}
}
```

---

## Methods Section Template

Use this template when describing MucOneSpan in your methods:

> MUC1 VNTR analysis was performed using MucOneSpan v[VERSION] (Popp, 2026; https://github.com/berntpopp/MucOneSpan), an open-source pipeline informed by some ideas from Vrbacka et al. (2025). PacBio HiFi CCS reads were mapped to a synthetic reference ladder (150 contigs, 1-150 repeat units) using minimap2. Allele lengths were determined using gap-based clustering with indel-valley splitting. Variants were called per allele using Clair3 with the PacBio HiFi model. Consensus sequences were classified into repeat units using the Vrbacka nomenclature with pre-computed mutation template matching.

---

## Acknowledging Tool Dependencies

MucOneSpan integrates several external tools. Please cite them as appropriate:

- **minimap2** -- Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics.* 2018;34(18):3094-3100.
- **Clair3** -- Zheng Z, et al. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. *Nat Comput Sci.* 2022;2:797-803.
- **samtools/bcftools** -- Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008.

---

## Contact

**Maintainer:** Bernt Popp

**GitHub:** [@berntpopp](https://github.com/berntpopp)

**Repository:** [berntpopp/MucOneSpan](https://github.com/berntpopp/MucOneSpan)

**Issues:** [GitHub Issue Tracker](https://github.com/berntpopp/MucOneSpan/issues)
