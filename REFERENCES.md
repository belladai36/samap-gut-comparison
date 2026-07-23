# Citation and data attribution

## Cite this repository

If you use the code, figures, or derived tables in this repository, cite the
repository using [`CITATION.cff`](CITATION.cff):

> Dai, J. (2026). *Pristina gut cross-species SAMap project*.
> https://github.com/belladai36/samap-gut-comparison

This repository is distributed under the proprietary terms in
[`LICENSE`](LICENSE). Permission is required before reuse, modification, or
redistribution.

## Upstream datasets

### Pristina leidyi

- Álvarez-Campos, P. *et al.* (2024). Annelid adult cell type diversity and
  their pluripotent cellular origins. *Nature Communications*, 15, 3194.
  https://doi.org/10.1038/s41467-024-47401-6
- Gene Expression Omnibus accession GSE230505:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230505
- Source atlas repository: https://github.com/scbe-lab/pristina-cell-type-atlas

### Human intestine

- Elmentaite, R. *et al.* (2021). Cells of the human intestinal tract mapped
  across space and time. *Nature*, 597, 250–255.
  https://doi.org/10.1038/s41586-021-03852-1
- Gut Cell Atlas portal: https://www.gutcellatlas.org/
- The analysis used the portal's epithelial raw-count AnnData download
  (`epi_raw_counts02_v2.h5ad`) and selected healthy adult epithelial/gut cells
  during preparation.

### Zebrafish intestine

- Jones, L. O. *et al.* (2023). Single-cell resolution of the adult zebrafish
  intestine under conventional conditions and in response to an acute *Vibrio
  cholerae* infection. *Cell Reports*, 42, 113407.
  https://doi.org/10.1016/j.celrep.2023.113407
- Gene Expression Omnibus accession GSE230044:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230044
- Broad Single Cell Portal study SCP2141:
  https://singlecell.broadinstitute.org/single_cell/study/SCP2141/adult-zebrafish-intestine
- Park, J. *et al.* (2019). Lysosome-rich enterocytes mediate protein
  absorption in the vertebrate gut. *Developmental Cell*, 51, 7–20.e6.
  https://doi.org/10.1016/j.devcel.2019.08.001

## Method and software

- Tarashansky, A. J. *et al.* (2021). Mapping single-cell atlases throughout
  Metazoa unravels cell type evolution. *eLife*, 10, e66747.
  https://doi.org/10.7554/eLife.66747
- SAMap source repository: https://github.com/atarashansky/SAMap
- NCBI BLAST+ command-line suite: Camacho, C. *et al.* (2009). BLAST+:
  architecture and applications. *BMC Bioinformatics*, 10, 421.
  https://doi.org/10.1186/1471-2105-10-421
- Scanpy: Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY:
  large-scale single-cell gene expression data analysis. *Genome Biology*, 19,
  15. https://doi.org/10.1186/s13059-017-1382-0
- AnnData: Virshup, I. *et al.* (2024). anndata: Access and store annotated
  data matrices. *Journal of Open Source Software*, 9, 4371.
  https://doi.org/10.21105/joss.04371
- Ensembl sequence resources: Yates, A. D. *et al.* (2020). Ensembl 2020.
  *Nucleic Acids Research*, 48, D682–D688.
  https://doi.org/10.1093/nar/gkz966

Human protein sequences were obtained from Ensembl release 98 and zebrafish
protein sequences from Ensembl release 96. Pristina transcript sequences and
annotations came from the Pristina source study and atlas repository above.

## Access and interpretation note

Dataset and portal links above were recorded on 2026-07-22. The public
repository contains derived tables and figures only; it does not redistribute
the upstream raw sequencing datasets. SAMap scores are cross-species
similarities, not probabilities or proof of one-to-one cell-type identity.
