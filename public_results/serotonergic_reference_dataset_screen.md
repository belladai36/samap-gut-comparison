# Serotonergic reference dataset screen

This screen evaluates candidate positive-reference datasets for future SAMap
score calibration. A dataset is considered immediately usable only when the
expression matrix, gene identifiers, per-cell biological labels, and sample
metadata can all be reproduced from public files.

## Human: GSE76381

- **Decision:** usable after deterministic gene-symbol-to-Ensembl translation.
- **Population:** the published `Sert` class from human fetal ventral midbrain.
- **Public data:** raw molecule-count matrix with `Cell_ID`, `Cell_type`, and
  developmental `Timepoint` metadata in the same CEF file.
- **Project conversion:** `prepare_human_serotonergic_reference.py` retains
  unique exact symbol matches to the Ensembl vocabulary already used by this
  project and stores the original counts in `X` and `layers["counts"]`.
- **Observed quality:** 1,977 cells and 16,105 matched Ensembl genes (82.5% gene
  coverage), with no duplicate identifiers, missing required metadata,
  negative counts, or non-integer counts.
- **Positive-class support:** 14 `hSert` cells from three embryo identifiers at
  weeks 9, 10, and 11. TPH2, SLC6A4, and FEV were detected in 64.3%, 71.4%,
  and 71.4% of `hSert` cells, compared with 0.6%, 2.1%, and 2.3% of other
  cells, respectively.
- **Limitation:** this is a developmental, plate-based dataset rather than an
  adult 10x atlas, and the positive class contains only 14 cells. It is suitable
  as an independent biological positive control, but is too small to serve as
  the sole training set for probability calibration. Batch, stage, and platform
  effects must be represented in sensitivity analyses.

## Zebrafish: GSE226494

- **Decision:** expression counts pass the initial format check; reference
  labels are pending.
- **Population:** published intraspinal serotonergic neurons (ISNs) from the
  `sox1a:eGFP` lineage at 1, 2, 3, and 5 dpf.
- **Public data:** four valid 10x filtered feature-barcode matrices.
- **Annotation evidence:** the authors' public code defines ISNs using their
  integrated Seurat clusters and reports markers including `tph2`, `slc6a4a`,
  `fev`, `lmx1bb`, `gata3`, and `gata2a`.
- **Blocking issue:** neither the GEO archive nor the authors' repository
  includes the final Seurat object or a barcode-to-cell-type table. Recreating
  clusters from raw matrices would be a reanalysis, not recovery of published
  labels, and therefore must not be treated as ground truth without an explicit
  reproducibility check or author-provided metadata.

## Next decision gate

1. Run the human conversion and quantify `Sert` representation, gene coverage,
   library coverage, sparsity, and marker detection.
2. Request or locate the zebrafish barcode-to-label metadata. If unavailable,
   rerun the authors' workflow and label the result as a reconstructed
   annotation with sensitivity analyses across clustering resolutions.
3. Add matched negative neuronal populations and reserve independent examples
   for probability-calibration evaluation.

## Sources

- Human GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76381>
- Human study: <https://pmc.ncbi.nlm.nih.gov/articles/PMC5055122/>
- Zebrafish GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226494>
- Zebrafish study: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10387610/>
- Zebrafish analysis code: <https://github.com/Fushun-Chen/Sox1a_scRNA>
