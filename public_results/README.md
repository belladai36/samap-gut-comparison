# Public SAMap result summary

This folder contains compact, shareable outputs from the completed three-species SAMap analysis (Pristina, human, zebrafish). It excludes raw sequencing data, `.h5ad` expression objects, BLAST databases, and the serialized SAMap object.

Use of this material requires written permission under the repository
[`LICENSE`](../LICENSE). Cite this project through [`CITATION.cff`](../CITATION.cff)
and credit all upstream datasets and tools listed in [`REFERENCES.md`](../REFERENCES.md).

## Files

- `pristina_cross_species_summary.csv`: best and second-best cell-type matches, scores, and score margins.
- `pristina_marker_homologs_top3.csv`: Pristina gut marker genes with the top three BLAST-supported human and zebrafish homolog candidates.
- `best_match_scores.png`: best match score for each Pristina gut cell type against each comparison species.
- `match_specificity_margins.png`: difference between the best and second-best match for each cell type and species.
- `priority_marker_validation.csv`: Ensembl-annotated top marker candidates for
  the three priority correspondences.
- `priority_candidate_assessment.md`: conservative interpretation of the
  priority marker evidence and the next validation experiments.
- `cross_species_module_enrichment.csv`: normalized, translated marker-module
  enrichment statistics for four pre-specified candidate correspondences.
- `cross_species_module_enrichment.png`: visualization of the cell-level
  permutation z-scores.
- `cross_species_module_enrichment_methods.md`: method and inferential
  limitations for the permutation analysis.
- `translated_reference_marker_modules.csv`: reference markers and their best
  BLAST-supported Pristina targets used in the module tests.
- `lre_validation_panel.csv` and `lre_validation_panel.md`: a five-gene,
  expression- and homology-screened panel for spatial testing of the LRE-like
  hypothesis.
- `lre_validation_panel.png`: expression enrichment and detection-specificity
  visualization for the five primary panel genes.
- `lre_mechanistic_context_checks.csv`: canonical zebrafish LRE machinery genes
  that are not primary positive Pristina probes under the current data.
- `lre_panel_cluster_specificity.csv` and `lre_panel_target_summary.csv`:
  cluster-by-gene expression and target-population summary for the five primary
  LRE-like candidates.
- `lre_panel_cluster_specificity.png` and `lre_panel_umap_expression.png`:
  the corresponding dot plot and computational expression display. A UMAP is
  not spatial tissue evidence.

## Interpretation

SAMap scores are cross-species similarity scores, not probabilities and not proof of one-to-one cell-type identity. The `specificity_class` column is a descriptive ranking heuristic based only on the top-score margin: higher specificity (>=0.20), moderate specificity (0.05–<0.20), or ambiguous (<0.05). Candidate matches should be validated with marker genes and biological context.
