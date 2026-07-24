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

## Interpretation

SAMap scores are cross-species similarity scores, not probabilities and not proof of one-to-one cell-type identity. The `specificity_class` column is a descriptive ranking heuristic based only on the top-score margin: higher specificity (>=0.20), moderate specificity (0.05–<0.20), or ambiguous (<0.05). Candidate matches should be validated with marker genes and biological context.
