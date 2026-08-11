# Professor-suggested reference cell types

## Outcome

- Pristina clusters 5 (`muscle 1`) and 6 (`muscle 2`) are supported striated-muscle reference candidates. Their published marker lists are dominated by myosin, a striated-muscle regulatory light chain, troponin, tropomyosin, and related contractile genes.
- Cluster 3 (`neurons 1`) contains the leading serotonergic candidate, but the evidence identifies a rare subpopulation rather than the entire cluster. All 46 neuronal cells with detectable `PrileiEVm011352t1` (annotated as TPH1) are in cluster 3. Of those, 9 also have detectable VMAT-like expression and 3 have detectable DDC-like expression.

## Important limitation

The current human dataset contains intestinal epithelial populations, and the current zebrafish dataset contains gut epithelial, immune, stromal, and pancreatic populations. Neither comparison dataset includes an annotated striated-muscle or neuronal population. Therefore, the existing gut-only SAMap run cannot provide a valid positive cross-species score for these suggested references.

## Next analysis

1. Subcluster the 5,138 cells in Pristina `neurons 1` and test whether TPH1-, VMAT-, and DDC-positive cells form a reproducible group across experimental libraries.
2. Acquire human and zebrafish reference datasets containing striated muscle and serotonergic neurons, with donor/sample metadata retained.
3. Prepare balanced `.h5ad` inputs and rerun SAMap on the expanded reference panel.
4. Use the muscle pairs as high-confidence positives. Keep the serotonin pair provisional until cluster specificity and marker coexpression are independently supported.
5. Construct negative pairs from biologically distinct types only after defining the positive reference populations.

## Sources

- Pristina atlas paper: https://doi.org/10.1038/s41467-024-47401-6
- GEO expression data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230505
- Supplementary Data 1: transcript annotation table.
- Supplementary Data 3: cluster differential-expression marker tables.
- Supplementary Data 4: top pseudobulk marker tables.

SAMap values are similarity scores, not calibrated probabilities. These candidate labels describe evidence available at this stage and should not be treated as confirmed evolutionary identity.
