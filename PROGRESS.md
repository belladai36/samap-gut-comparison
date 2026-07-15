# Project progress

Last updated: 2026-07-15

## Completed

- Created a project-local Python 3.12 environment.
- Installed and verified SAMap 3.0.1, Scanpy, AnnData, Snakemake, and NCBI BLAST.
- Defined a three-species configuration template.
- Implemented fail-fast validation for expression matrices, metadata, protein
  FASTAs, identifier overlap, checksums, and homology maps.
- Implemented an official SAMap baseline stage using the installed 3.0.1 API.
- Implemented a separate LocalSAMap stage with blocked permutations,
  Benjamini-Hochberg q-values, specificity margins, and mapping entropy.
- Verified workflow construction with a complete Snakemake dry-run.
- Passed a synthetic LocalSAMap graph-aggregation smoke test.

## Current blocker

The analysis cannot run until the configured expression datasets, matching
protein FASTAs, and SAMap BLAST homology maps are available locally. These
inputs are intentionally excluded from version control.

## Next milestone

Prepare and validate the focal-species expression matrix and matching protein
FASTA, then add the two comparison datasets and generate reciprocal BLAST maps.

