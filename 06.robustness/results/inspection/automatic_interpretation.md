
### Descriptive result summary

- The analysis contains **369 proteins** evaluated over **20 repeated 80/20 splits**.
- At Pearson ≥ **0.50**, required in at least **80%** of splits, **2 (0.5%)** proteins were classified as proxy-sufficient and **360 (97.6%)** as robustly linear-improvable.
- **7 (1.9%)** proteins were unstable or borderline, and **0 (0.0%)** were unresolved under the saved rule.
- The median of the per-protein median Pearson values was **0.312** for the proxy and **0.702** for ElasticNet.
- The median repeated-split ElasticNet-minus-proxy gain was **0.387**. The descriptive gain interval lower bound exceeded **0.05** for **369** proteins.

### Recommended wording for the manuscript

> Across 20 repeated 80/20 partitions of CCLE, we quantified protein-level stability using the median and interquartile range of held-out Pearson correlations and the fraction of partitions exceeding a prespecified correlation threshold. At Pearson ≥ 0.50 in at least 80% of partitions, 2 proteins were consistently captured by a fitted same-gene RNA proxy, whereas 360 additional proteins were consistently captured by ElasticNet and showed a descriptive median gain of at least 0.05 whose split-level bootstrap lower bound also exceeded that value.

### Required caveat

> Repeated holdout test sets overlap, and the bootstrap intervals summarize the empirical distribution across repeated partitions rather than independent experimental replicates. We therefore interpret this analysis as a stability assessment and use full-sample out-of-fold and external-cohort evaluations for primary inferential and transfer claims.
