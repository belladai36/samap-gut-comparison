# Pristina LRE-like spatial validation panel

This panel prioritizes genes for testing the hypothesis that Pristina anterior/mid-intestine 2 contains an LRE-like program.
It is not a definitive cell-type annotation.

## Selection logic

- Target population: `anterior/mid-intestine 2`.
- Primary positive panel: lysosomal/protein-digestion LRE markers with high BLAST alignment support, positive Pristina target-versus-other expression, and >0.05 detection-fraction difference.
- Context checks: zebrafish `cubn`, `dab2`, and `amn` are retained separately because they form LRE-associated receptor-mediated endocytic machinery, but they are not primary positive probes unless they are enriched in Pristina.
- Pristina expression is library-size normalized to 10,000 counts per cell and log1p transformed before scoring.
- Biological rationale source: Park et al. 2019, Developmental Cell, doi:10.1016/j.devcel.2019.08.001.

## Recommended experiment

Use a multiplex panel of several primary positive genes. Confirm that signals co-localize in anterior/mid-intestine 2, then compare the localization with established gut-region markers. The canonical zebrafish machinery genes should be interpreted as separate context checks, not assumed to be positive Pristina markers. A negative or discordant result should be interpreted as evidence against, not merely absence of proof for, the LRE-like hypothesis.

## Panel

| Rank | Zebrafish symbol | Pristina gene | Class | Homology quality | Identity (%) | Alignment length | Pristina log2 FC | Target detection | Detection difference |
|---:|---|---|---|---|---:|---:|---:|---:|---:|
| 1 | fuca1.1 | PrileiEVm011523t1 | lysosomal/protein-digestion support | high | 57.143 | 455 | 3.712 | 0.580 | 0.554 |
| 2 | hexb | PrileiEVm008813t1 | lysosomal/protein-digestion support | high | 56.897 | 522 | 3.667 | 0.547 | 0.522 |
| 3 | ctsl.1 | PrileiEVm014063t1 | lysosomal/protein-digestion support | high | 64.286 | 294 | 2.122 | 0.358 | 0.289 |
| 4 | ctsbb | PrileiEVm014379t1 | lysosomal/protein-digestion support | high | 62.379 | 311 | 2.660 | 0.187 | 0.175 |
| 5 | naga | PrileiEVm012889t1 | lysosomal/protein-digestion support | high | 52.507 | 379 | 2.777 | 0.237 | 0.221 |

## Caveats

- BLAST-supported homology is not proof of one-to-one orthology or conserved function.
- The expression statistics are single-cell descriptive quantities, not replicate-level inference.
- `cubn`, `dab2`, and `amn` are reported separately in `lre_mechanistic_context_checks.csv`; they are not used as primary positive probes because their current Pristina expression does not support that use.
