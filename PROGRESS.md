# Living project report: Pristina gut cross-species comparison

**Last updated:** 2026-07-28  
**Project status:** baseline SAMap analysis and exploratory module validation complete

## Aim

Compare annotated Pristina *leidyi* gut populations with human and zebrafish
intestinal cell atlases. The analysis uses SAMap for cross-species cell-state
similarity and BLAST homology maps to connect non-identical gene repertoires.

## Completed work

1. Created a project-local Python 3.12 analysis environment containing SAMap,
   Scanpy, AnnData, and NCBI BLAST+.
2. Prepared three input expression objects:
   - Pristina gut: 9,253 cells × 27,758 genes.
   - Human gut (balanced input): 8,732 cells × 33,538 genes.
   - Zebrafish gut: 18,358 cells × 25,107 genes.
3. Generated and quality-checked reciprocal BLAST tables for Pristina–human,
   Pristina–zebrafish, and human–zebrafish gene homology.
4. Ran the official SAMap workflow for all three species and saved the joint
   mapping object locally. Large raw/intermediate files are deliberately not
   version-controlled.
5. Exported compact mapping summaries and marker-homology tables to
   [`public_results/`](public_results/).
6. Performed a pre-specified, exploratory module-enrichment validation using
   normalized expression, translated reference markers, 2,000 label
   permutations, and Benjamini–Hochberg FDR correction.
7. Built a spatial-validation priority panel for the LRE-like hypothesis from
   the translated LRE module, Pristina expression specificity, and alignment
   quality. The canonical zebrafish `cubn`, `dab2`, and `amn` candidates were
   retained as context checks, not primary positive probes, because they are
   not enriched in the Pristina target population.
8. Confirmed the five primary LRE-like panel genes across every annotated
   Pristina gut population. `anterior/mid-intestine 2` has rank 1 mean
   normalized expression for all five genes (`fuca1.1`, `hexb`, `ctsl.1`,
   `ctsbb`, and `naga`). This is cluster-specificity support, not spatial
   validation.
9. Prepared primary and backup 400-nt probe regions for all five LRE-like
   targets. Each recommended region passed a conservative exact-20-mer
   transcriptome screen with zero strong non-self matches. These are
   platform-neutral target regions, not final assay oligos.

## Current findings

| Pristina population | Reference candidate | Supporting evidence | Interpretation |
|---|---|---|---|
| anterior/mid-intestine 2 | zebrafish LRE | SAMap score 0.762; module z = 47.3; FDR = 0.0005 | Strongest current candidate; consistent with a lysosome-rich-enterocyte-like program. |
| anterior/mid-intestine 2 | human Paneth | SAMap score 0.377; module z = 12.4; FDR = 0.0005 | A conserved antimicrobial/lysosomal program is plausible, but this does not establish Paneth identity. |
| anterior intestine | human enterocyte | SAMap score 0.445; module z = 6.7; FDR = 0.0005 | Supports an absorptive/enterocyte-like hypothesis. |
| stomach 2 | human colonocyte | SAMap score 0.593; module z = 5.4; FDR = 0.0005 | A candidate similarity; current top markers are not colonocyte-specific enough for a definitive label. |

SAMap scores are similarity scores, not probabilities or one-to-one cell-type
identity calls. The permutation tests are **cell-level exploratory tests**;
they do not turn cells into independent biological replicates or account for
all donor, batch, or evolutionary dependence. The zebrafish input has one
sample, so its inference must be especially cautious.

For the LRE-like hypothesis, the primary spatial panel comprises five
Pristina-enriched candidates connected to zebrafish lysosomal/protein-digestion
genes: `PrileiEVm011523t1` (`fuca1.1`), `PrileiEVm008813t1` (`hexb`),
`PrileiEVm014063t1` (`ctsl.1`), `PrileiEVm014379t1` (`ctsbb`), and
`PrileiEVm012889t1` (`naga`). This supports testing a **lysosomal program**;
it does not prove complete conservation of canonical zebrafish LRE machinery.

## Main outputs

- [Cross-species best-match summary](public_results/pristina_cross_species_summary.csv)
- [Priority marker assessment](public_results/priority_candidate_assessment.md)
- [Module-enrichment statistics](public_results/cross_species_module_enrichment.csv)
- [Module-enrichment methods and limitations](public_results/cross_species_module_enrichment_methods.md)
- [Marker modules translated through BLAST](public_results/translated_reference_marker_modules.csv)
- [LRE-like spatial validation panel](public_results/lre_validation_panel.md)
- [LRE panel cluster-specificity table](public_results/lre_panel_cluster_specificity.csv)
- [LRE panel target summary](public_results/lre_panel_target_summary.csv)
- [LRE probe-region design and limitations](public_results/lre_probe_region_design.md)
- [Recommended probe regions](public_results/lre_probe_region_recommendations.csv)
- [LRE mechanistic context checks](public_results/lre_mechanistic_context_checks.csv)
- [Public result figures and tables](public_results/README.md)

## Reproducible notebooks

- [01 — SAMap mapping summary](notebooks/01_samap_mapping_summary.ipynb)
- [02 — Cross-species module validation](notebooks/02_cross_species_module_validation.ipynb)
- [03 — Priority marker review](notebooks/03_priority_marker_review.ipynb)
- [04 — LRE-like spatial validation panel](notebooks/04_lre_validation_panel.ipynb)

Each notebook uses the shareable tables in `public_results/`; it therefore does
not require the large `.h5ad` objects or serialized SAMap model. Run from the
repository root with the project environment:

```bash
"$PWD/.conda/envs/samap-gut/bin/python" -m jupyter nbconvert \
  --execute --to notebook --inplace notebooks/01_samap_mapping_summary.ipynb
```

## Next steps

1. Select an assay platform (HCR RNA-FISH, smFISH, or ISH/RNAscope) and use the
   recommended target regions to create final assay-specific oligos.
2. Validate the leading candidates using spatial expression or in situ methods;
   require co-localization in anatomical anterior/mid intestine rather than
   interpreting UMAP position as tissue location.
3. When independent biological replicates are available, repeat the tests with
   donor/animal-aware pseudobulk rather than cell-level label permutations.
4. Compare original SAMap outputs with the separate local-pattern method only
   after defining the comparison endpoints and validation criteria in advance.

## How to keep this report current

When a new analysis milestone is completed:

1. Update the **Last updated** date and append the result under **Completed
   work** or **Current findings**.
2. Add any new small, shareable result table or figure to `public_results/` and
   link it under **Main outputs**.
3. Add or update a notebook if the result needs code, a figure, and a written
   interpretation.
4. Record methodological changes and limitations rather than replacing older
   conclusions silently.

## Data handling

Raw expression matrices, BLAST databases, full mapping objects, and other large
intermediates remain local and excluded from Git. Dataset and software sources
are listed in [REFERENCES.md](REFERENCES.md); usage conditions are in
[LICENSE](LICENSE).
