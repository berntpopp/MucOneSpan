# Literature check for the two length failures

GeneFoundry MCP/PubTator retrieved2026-09-14. Raw search and full METHODS passages
are in PV/phasing_debug/genefoundry_length_model_search.json and
genefoundry_nanorepeat_methods.json. This is method precedent, not MUC1 validation.

**Recommended citation:** Fang L et al. Haplotyping SNPs for allele-specific gene editing of the expanded huntingtin allele using long-read sequencing. HGG Adv. 2023. PMID:36262216. doi:10.1016/j.xhgg.2022.100146.

[Primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9574884/).
Passages PMID:36262216:methods:15 and :16 describe NanoRepeat: per-read repeat
estimates by reference-fit refinement, followed by a one/two-component Gaussian
mixture. The authors found model-selection overfitting on integer-rounded counts
and introduced random jitter for BIC selection. They also removed global three-SD
outliers and excluded samples with one component or fewer than50reads per allele.

These are unsuitable defaults for this task without validation: global outlier
removal can delete the real long minority allele, a50-read cutoff rejects the
explicit low-coverage objective, and excluding one component cannot establish
sequence homozygosity. Random jitter would require a configured recorded seed and
would be avoided if continuous base-span likelihood is used. Continuous spans
still have uncertain boundaries, read error, and length-dependent variance;
BIC is not a calibrated genotype confidence or guarantee against false splitting.

Together with the earlier HMMSTR/Straglr evidence in
genefoundry-pubtator-research.md, this motivates testing a bounded likelihood
alternative against all exposed controls. It does not justify promoting a GMM
because it repairs two selected examples. The existing no-regression, exact count,
assignment and downstream mutation/reconstruction gates remain unchanged.
