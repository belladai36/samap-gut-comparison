# Priority candidate assessment

This is a focused marker-based check of the three highest-priority SAMap
correspondence candidates. It uses the top 15 Pristina marker genes per cell
type when a top-ranked BLAST candidate could be annotated by Ensembl. See
`priority_marker_validation.csv` for the individual rows.

## 1. Anterior/mid-intestine 2 → zebrafish LRE

This is the most strongly supported candidate at present. The SAMap result is
high and specific (score 0.762; best-minus-second margin 0.731). Its top
annotated marker candidates include lysosomal/catabolic genes: `CTSL`, `FUCA2`,
`HEXB`, `MAN2B1`, and `GALC`, together with repeated `MRC1` and `CPAMD8`
candidates. This is compatible with a lysosome-rich program.

In zebrafish, LRE means **lysosome-rich enterocyte**: an absorptive epithelial
cell type with strong endocytic and lysosomal activity. The LRE interpretation
is therefore a biologically coherent hypothesis, not a transferred final name.
It should next be tested with a broader LRE marker set, especially `LAMP2`,
`CUBN`, `AMN`, and `DAB2`, plus spatial/anatomical localization.

The parallel human Paneth match (score 0.377; margin 0.301) remains a separate
candidate. The current top-marker subset does not independently establish a
canonical Paneth program, so do not label this Pristina population as Paneth
without additional evidence.

## 2. Stomach 2 → human colonocyte

The SAMap evidence is strong and specific (score 0.593; margin 0.564), making
this a high-priority correspondence to validate. The top annotatable marker
candidates include `ABCB1`, but many of the remaining markers are mitochondrial
or broadly expressed energy/metabolism genes (`MT-CO1`, `SLC25A3`, `ATP5F1A`,
`ATP5F1B`, `ACTB`). This compact marker subset does **not** yet supply a
colonocyte-specific molecular signature.

Next, compare a larger marker set against human colonocyte programs and assess
whether genes related to epithelial transport, ion handling, and absorptive
metabolism are consistently enriched. Keep the label as “human-colonocyte-like
candidate” until this check is complete.

## 3. Anterior intestine → human enterocyte

This correspondence is reasonably separated in the human comparison (score
0.445; margin 0.207). Annotated candidates include `FAR1`, `PLB1`, `ACSBG1`,
and `SLC6A9`, consistent with a broad metabolic/transport-oriented epithelial
program, but not uniquely diagnostic of enterocytes in this small panel.

Next, test a full absorptive-enterocyte gene set rather than relying on
individual BLAST hits. A stable conclusion should agree across: (1) SAMap
score/margin, (2) a multi-gene marker/module score, and (3) anatomical position
or spatial validation.

## Limits and next experiment

Each listed homolog is the top candidate from the project's BLAST-derived map;
it is not automatically a confirmed one-to-one ortholog. BLAST evidence,
single-cell marker enrichment, and SAMap mappings should be interpreted jointly.
The recommended next computational test is a pre-specified module score for
LRE, enterocyte, colonocyte, and Paneth marker sets, followed by a permutation
test at the cell-type/pseudobulk level. Experimental localization of selected
Pristina markers is needed before final anatomical label transfer.

## Sources

- SAMap mapping and score interpretation: Tarashansky *et al.* (2021),
  https://doi.org/10.7554/eLife.66747
- LRE definition and lysosomal/endocytic program: Park *et al.* (2019),
  https://doi.org/10.1016/j.devcel.2019.08.001
- Dataset and software attributions: [`../REFERENCES.md`](../REFERENCES.md)
