# Cross-species module-enrichment test

This exploratory analysis tests four candidate correspondences selected from the
SAMap mapping results. For each human or zebrafish reference cell type, we
normalized every count matrix to 10,000 counts per cell and log1p-transformed
it. We then derived up-regulated markers with Scanpy's `t-test_overestim_var` comparison
against all other reference cell types (adjusted p < 0.05; positive log fold
change; at most 50 genes). Each source gene was translated
to the best Pristina BLAST hit by minimum E-value. The Pristina score is the
mean gene-wise z-scored expression over the translated module.

For a pre-specified Pristina population of size n, the statistic is its mean
module score. The one-sided p-value uses 2,000 random label-set
permutations of size n:

`p = (1 + #{T_perm >= T_obs}) / (2000 + 1)`.

BH FDR is applied across the four planned tests. This is a *cell-level,
exploratory* test: its label-exchangeability null does not make cells biological
replicates, and it does not control for donor, batch, or evolutionary dependence.
In particular, the zebrafish dataset has one sample. Candidate identities still
require independent biological validation, preferably donor/animal-aware
pseudobulk and spatial or in situ evidence.
