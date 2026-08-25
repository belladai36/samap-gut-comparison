# Pristina serotonergic-neuron subclustering

## Outcome

At Leiden resolution 0.6, cluster `6` is the TPH1-enriched candidate. It contains 226 of 5,138 `neurons 1` cells and captures 30 of 46 TPH1-like-positive cells. TPH1-like detection is 13.27% in the candidate and 0.33% in the remaining cells (one-sided Fisher exact test, BH-adjusted P = 3.65e-30). A 10,000-permutation test correcting for selection of the most enriched cluster gives empirical P = 9.999e-05 (0 shuffled datasets were as or more extreme). TPH1-positive candidate cells occur in 7 of 9 libraries from 2 of 3 experiments.

**Assessment: supported but still provisional.** This is evidence for a rare serotonergic-like subpopulation, not proof of a conserved cell-type identity. Sparse transcript detection, dropout, and clustering sensitivity remain important limitations.

## Validation criteria

- Marker enrichment: one-sided Fisher exact tests with Benjamini-Hochberg correction.
- Cluster selection: a fixed-resolution permutation test compares the best observed cluster with the best cluster under shuffled TPH labels.
- Resolution stability: TPH1-enriched communities were compared from Leiden resolutions 0.2 to 0.8 using cell-set Jaccard similarity.
- Reproducibility: candidate and TPH1-positive counts were tabulated separately for all nine libraries and three experiments.
- Differential expression: candidate versus all other `neurons 1` cells was tested using Scanpy's Wilcoxon implementation.

## Files

- `serotonergic_subcluster_summary.csv`: headline candidate statistics.
- `serotonergic_marker_enrichment.csv`: marker detection and Fisher tests.
- `serotonergic_resolution_stability.csv`: sensitivity to Leiden resolution.
- `serotonergic_cluster_selection_permutation.csv`: post-selection permutation test.
- `serotonergic_library_support.csv`: experiment/library-level support.
- `serotonergic_candidate_markers.csv`: top differential markers.
- `serotonergic_neuron_subclustering_umap.png`: clusters and marker expression.
- `serotonergic_marker_enrichment.png`: candidate-versus-other detection fractions.

## Next decision

Use this candidate as a provisional Pristina query population only after reviewing its differential markers and confirming it is not driven by one library. Cross-species probability calibration still requires vertebrate serotonergic reference populations and curated negative controls.
