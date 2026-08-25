"""Subcluster Pristina ``neurons 1`` and evaluate serotonergic evidence.

The analysis is deliberately conservative: a Leiden community is selected by
enrichment for the TPH1-like transcript, but it is not called a confirmed cell
type.  Coexpression, resolution stability, and support across libraries are
reported separately so that the evidence can be audited.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/raw/pristina/pristina_atlas_final.h5ad"
RESULTS = ROOT / "results/serotonergic_neuron_subclustering"
PUBLIC = ROOT / "public_results"

SEED = 230505
TARGET_PARENT_CLUSTER = "3"
TARGET_PARENT_LABEL = "neurons 1"
SELECTED_RESOLUTION = 0.6
RESOLUTIONS = (0.2, 0.4, 0.6, 0.8)

MARKERS = {
    "TPH1-like": "PrileiEVm011352t1",
    "DDC-like": "PrileiEVm010940t1",
    "VMAT-like 1": "PrileiEVm008386t1",
    "VMAT-like 2": "PrileiEVm008853t1",
    "DAT-like supporting": "PrileiEVm006035t1",
}


def _dense(vector) -> np.ndarray:
    if sparse.issparse(vector):
        return vector.toarray().ravel()
    return np.asarray(vector).ravel()


def _bh(values: list[float]) -> np.ndarray:
    return multipletests(values, method="fdr_bh")[1]


def prepare_neurons() -> sc.AnnData:
    atlas = sc.read_h5ad(SOURCE)
    parent = atlas.obs["leiden_1.5"].astype(str) == TARGET_PARENT_CLUSTER
    neurons = atlas[parent].copy()
    del atlas

    # Reconstruct preprocessing from integer-like raw counts in the retained
    # feature set rather than reclustering the atlas's scaled expression.
    neurons.X = neurons.layers["counts"].copy()
    neurons.layers["counts"] = neurons.X.copy()
    sc.pp.filter_genes(neurons, min_cells=10)
    sc.pp.normalize_total(neurons, target_sum=1e4)
    sc.pp.log1p(neurons)
    neurons.raw = neurons

    sc.pp.highly_variable_genes(neurons, n_top_genes=3000, flavor="seurat")
    missing_markers = sorted(set(MARKERS.values()) - set(neurons.var_names))
    if missing_markers:
        raise ValueError(f"Candidate genes absent after filtering: {missing_markers}")
    # Do not force the hypothesis markers into the clustering feature set.
    # Marker enrichment is evaluated only after unsupervised clustering, which
    # prevents TPH1/DDC/VMAT expression from defining the candidate by design.
    neurons = neurons[:, neurons.var["highly_variable"]].copy()
    sc.pp.scale(neurons, max_value=10)
    sc.tl.pca(neurons, n_comps=50, random_state=SEED)
    sc.pp.neighbors(neurons, n_neighbors=20, n_pcs=40, random_state=SEED)
    sc.tl.umap(neurons, random_state=SEED)
    for resolution in RESOLUTIONS:
        sc.tl.leiden(
            neurons,
            resolution=resolution,
            key_added=f"serotonin_leiden_{resolution:g}",
            random_state=SEED,
        )
    return neurons


def add_marker_detection(neurons: sc.AnnData) -> None:
    for label, gene in MARKERS.items():
        values = _dense(neurons.raw[:, gene].X)
        neurons.obs[f"detected_{label}"] = values > 0
        neurons.obs[f"expression_{label}"] = values
    neurons.obs["detected_any_VMAT"] = (
        neurons.obs["detected_VMAT-like 1"] | neurons.obs["detected_VMAT-like 2"]
    )
    neurons.obs["detected_TPH1_support"] = (
        neurons.obs["detected_TPH1-like"]
        & (
            neurons.obs["detected_DDC-like"]
            | neurons.obs["detected_any_VMAT"]
            | neurons.obs["detected_DAT-like supporting"]
        )
    )


def choose_tph_cluster(neurons: sc.AnnData, key: str) -> str:
    grouped = neurons.obs.groupby(key, observed=True)["detected_TPH1-like"].agg(["sum", "count"])
    grouped["fraction"] = grouped["sum"] / grouped["count"]
    # Detection fraction is primary; detected-cell count breaks ties and
    # prevents arbitrary selection among empty clusters.
    return str(grouped.sort_values(["fraction", "sum"], ascending=False).index[0])


def marker_enrichment(neurons: sc.AnnData, cluster_key: str, candidate: str) -> pd.DataFrame:
    inside = neurons.obs[cluster_key].astype(str) == candidate
    rows = []
    for label in [*MARKERS, "any_VMAT", "TPH1_support"]:
        detected = neurons.obs[f"detected_{label}"].to_numpy(dtype=bool)
        a = int(np.sum(detected & inside))
        b = int(np.sum(~detected & inside))
        c = int(np.sum(detected & ~inside))
        d = int(np.sum(~detected & ~inside))
        odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "marker": label,
                "candidate_positive": a,
                "candidate_total": a + b,
                "candidate_detection_fraction": a / (a + b),
                "other_positive": c,
                "other_total": c + d,
                "other_detection_fraction": c / (c + d),
                "odds_ratio": odds,
                "fisher_pvalue": pvalue,
            }
        )
    result = pd.DataFrame(rows)
    result["fdr_bh"] = _bh(result["fisher_pvalue"].tolist())
    return result


def resolution_stability(neurons: sc.AnnData, selected_mask: np.ndarray) -> pd.DataFrame:
    tph = neurons.obs["detected_TPH1-like"].to_numpy(dtype=bool)
    rows = []
    for resolution in RESOLUTIONS:
        key = f"serotonin_leiden_{resolution:g}"
        candidate = choose_tph_cluster(neurons, key)
        mask = neurons.obs[key].astype(str).to_numpy() == candidate
        intersection = int(np.sum(mask & selected_mask))
        union = int(np.sum(mask | selected_mask))
        rows.append(
            {
                "resolution": resolution,
                "cluster_key": key,
                "tph_enriched_cluster": candidate,
                "cluster_cells": int(mask.sum()),
                "tph_positive_cells": int(np.sum(tph & mask)),
                "tph_detection_fraction": float(np.mean(tph[mask])),
                "tph_recall": float(np.sum(tph & mask) / max(1, tph.sum())),
                "jaccard_vs_resolution_0.6": intersection / union,
            }
        )
    return pd.DataFrame(rows)


def cluster_selection_permutation(
    neurons: sc.AnnData, cluster_key: str, n_permutations: int = 10_000
) -> pd.DataFrame:
    """Correct the TPH enrichment test for selecting the best cluster.

    For each shuffled TPH label vector, the smallest one-sided hypergeometric
    P value across all clusters is recorded.  Comparing the observed minimum to
    this null distribution controls the fixed-resolution cluster search.
    """
    labels = neurons.obs[cluster_key].astype(str).to_numpy()
    tph = neurons.obs["detected_TPH1-like"].to_numpy(dtype=bool)
    clusters, codes = np.unique(labels, return_inverse=True)
    cluster_sizes = np.bincount(codes)
    total_cells = len(labels)
    total_tph = int(tph.sum())

    observed_counts = np.bincount(codes, weights=tph.astype(int)).astype(int)
    observed_pvalues = hypergeom.sf(observed_counts - 1, total_cells, total_tph, cluster_sizes)
    observed_minimum = float(np.min(observed_pvalues))

    rng = np.random.default_rng(SEED)
    tph_indices = np.arange(total_cells)
    null_minimum = np.empty(n_permutations, dtype=float)
    for iteration in range(n_permutations):
        selected = rng.choice(tph_indices, size=total_tph, replace=False)
        counts = np.bincount(codes[selected], minlength=len(clusters))
        pvalues = hypergeom.sf(counts - 1, total_cells, total_tph, cluster_sizes)
        null_minimum[iteration] = np.min(pvalues)
    exceedances = int(np.sum(null_minimum <= observed_minimum))
    empirical = (1 + exceedances) / (n_permutations + 1)
    return pd.DataFrame(
        [
            {
                "cluster_key": cluster_key,
                "clusters_searched": len(clusters),
                "n_permutations": n_permutations,
                "observed_minimum_hypergeometric_p": observed_minimum,
                "null_as_or_more_extreme": exceedances,
                "empirical_cluster_selection_p": empirical,
                "null_minimum_p_05": float(np.quantile(null_minimum, 0.05)),
                "null_minimum_p_50": float(np.quantile(null_minimum, 0.50)),
            }
        ]
    )


def library_support(neurons: sc.AnnData, selected_mask: np.ndarray) -> pd.DataFrame:
    frame = neurons.obs[["Experiment", "Library", "detected_TPH1-like", "detected_TPH1_support"]].copy()
    frame["candidate"] = selected_mask
    rows = []
    for (experiment, library), group in frame.groupby(["Experiment", "Library"], observed=True):
        candidate = group["candidate"].to_numpy(dtype=bool)
        tph = group["detected_TPH1-like"].to_numpy(dtype=bool)
        support = group["detected_TPH1_support"].to_numpy(dtype=bool)
        rows.append(
            {
                "experiment": str(experiment),
                "library": str(library),
                "neurons1_cells": len(group),
                "candidate_cells": int(candidate.sum()),
                "candidate_fraction": float(candidate.mean()),
                "tph_positive_cells": int(tph.sum()),
                "tph_positive_in_candidate": int(np.sum(tph & candidate)),
                "tph_supported_in_candidate": int(np.sum(support & candidate)),
            }
        )
    return pd.DataFrame(rows)


def differential_markers(neurons: sc.AnnData, cluster_key: str, candidate: str) -> pd.DataFrame:
    work = neurons.raw.to_adata()
    work.obs["candidate_serotonergic"] = np.where(
        neurons.obs[cluster_key].astype(str).to_numpy() == candidate, "candidate", "other_neurons1"
    )
    sc.tl.rank_genes_groups(
        work,
        groupby="candidate_serotonergic",
        groups=["candidate"],
        reference="other_neurons1",
        method="wilcoxon",
        pts=True,
    )
    markers = sc.get.rank_genes_groups_df(work, group="candidate").head(100)
    annotation_columns = [column for column in ["Description", "Preferred_name"] if column in work.var]
    if annotation_columns:
        annotations = work.var[annotation_columns].copy()
        annotations.index.name = "names"
        markers = markers.merge(annotations.reset_index(), on="names", how="left")
    return markers


def save_figures(neurons: sc.AnnData, cluster_key: str, candidate: str) -> None:
    PUBLIC.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = PUBLIC
    sc.pl.umap(
        neurons,
        color=[cluster_key, "expression_TPH1-like", "expression_DDC-like", "expression_VMAT-like 1"],
        cmap="viridis",
        ncols=2,
        show=False,
    )
    plt.savefig(PUBLIC / "serotonergic_neuron_subclustering_umap.png", dpi=220, bbox_inches="tight")
    plt.close("all")

    plot = neurons.obs.copy()
    plot["candidate_label"] = np.where(
        plot[cluster_key].astype(str) == candidate, "TPH1-enriched candidate", "other neurons 1"
    )
    order = ["other neurons 1", "TPH1-enriched candidate"]
    values = []
    for group in order:
        mask = plot["candidate_label"] == group
        values.append(
            [
                plot.loc[mask, "detected_TPH1-like"].mean(),
                plot.loc[mask, "detected_DDC-like"].mean(),
                plot.loc[mask, "detected_any_VMAT"].mean(),
                plot.loc[mask, "detected_TPH1_support"].mean(),
            ]
        )
    fig, ax = plt.subplots(figsize=(8, 4.8))
    x = np.arange(4)
    width = 0.36
    for index, group in enumerate(order):
        ax.bar(x + (index - 0.5) * width, values[index], width, label=group)
    ax.set_xticks(x, ["TPH1-like", "DDC-like", "any VMAT-like", "TPH1 + support"])
    ax.set_ylabel("Fraction of cells with detected expression")
    ax.set_title("Serotonergic marker detection in Pristina neurons 1")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(PUBLIC / "serotonergic_marker_enrichment.png", dpi=220)
    plt.close(fig)


def write_report(
    neurons: sc.AnnData,
    candidate: str,
    enrichment: pd.DataFrame,
    stability: pd.DataFrame,
    libraries: pd.DataFrame,
    selection_test: pd.DataFrame,
) -> None:
    key = f"serotonin_leiden_{SELECTED_RESOLUTION:g}"
    mask = neurons.obs[key].astype(str).to_numpy() == candidate
    tph = neurons.obs["detected_TPH1-like"].to_numpy(dtype=bool)
    tph_row = enrichment.loc[enrichment["marker"] == "TPH1-like"].iloc[0]
    n_supported_libraries = int((libraries["tph_positive_in_candidate"] > 0).sum())
    n_supported_experiments = int(
        libraries.loc[libraries["tph_positive_in_candidate"] > 0, "experiment"].nunique()
    )
    min_jaccard = float(stability["jaccard_vs_resolution_0.6"].min())
    selection_p = float(selection_test["empirical_cluster_selection_p"].iloc[0])
    assessment = "provisional"
    if n_supported_libraries >= 3 and tph_row["fdr_bh"] < 0.05 and min_jaccard >= 0.5:
        assessment = "supported but still provisional"
    report = f"""# Pristina serotonergic-neuron subclustering

## Outcome

At Leiden resolution {SELECTED_RESOLUTION:g}, cluster `{candidate}` is the TPH1-enriched candidate. It contains {mask.sum():,} of {neurons.n_obs:,} `neurons 1` cells and captures {np.sum(mask & tph)} of {tph.sum()} TPH1-like-positive cells. TPH1-like detection is {tph_row['candidate_detection_fraction']:.2%} in the candidate and {tph_row['other_detection_fraction']:.2%} in the remaining cells (one-sided Fisher exact test, BH-adjusted P = {tph_row['fdr_bh']:.3g}). A {int(selection_test['n_permutations'].iloc[0]):,}-permutation test correcting for selection of the most enriched cluster gives empirical P = {selection_p:.4g} ({int(selection_test['null_as_or_more_extreme'].iloc[0])} shuffled datasets were as or more extreme). TPH1-positive candidate cells occur in {n_supported_libraries} of {len(libraries)} libraries from {n_supported_experiments} of 3 experiments.

**Assessment: {assessment}.** This is evidence for a rare serotonergic-like subpopulation, not proof of a conserved cell-type identity. Sparse transcript detection, dropout, and clustering sensitivity remain important limitations.

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
"""
    (PUBLIC / "serotonergic_neuron_subclustering.md").write_text(report)


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)
    RESULTS.mkdir(parents=True, exist_ok=True)
    PUBLIC.mkdir(parents=True, exist_ok=True)

    neurons = prepare_neurons()
    add_marker_detection(neurons)
    key = f"serotonin_leiden_{SELECTED_RESOLUTION:g}"
    candidate = choose_tph_cluster(neurons, key)
    selected = neurons.obs[key].astype(str).to_numpy() == candidate

    enrichment = marker_enrichment(neurons, key, candidate)
    stability = resolution_stability(neurons, selected)
    selection_test = cluster_selection_permutation(neurons, key)
    libraries = library_support(neurons, selected)
    markers = differential_markers(neurons, key, candidate)
    summary = pd.DataFrame(
        [
            {
                "parent_population": TARGET_PARENT_LABEL,
                "selected_resolution": SELECTED_RESOLUTION,
                "candidate_cluster": candidate,
                "neurons1_cells": neurons.n_obs,
                "candidate_cells": int(selected.sum()),
                "tph_positive_total": int(neurons.obs["detected_TPH1-like"].sum()),
                "tph_positive_in_candidate": int(
                    np.sum(selected & neurons.obs["detected_TPH1-like"].to_numpy(dtype=bool))
                ),
                "libraries_with_tph_candidate_cells": int((libraries["tph_positive_in_candidate"] > 0).sum()),
                "experiments_with_tph_candidate_cells": int(
                    libraries.loc[libraries["tph_positive_in_candidate"] > 0, "experiment"].nunique()
                ),
                "cluster_selection_empirical_p": float(
                    selection_test["empirical_cluster_selection_p"].iloc[0]
                ),
            }
        ]
    )

    summary.to_csv(PUBLIC / "serotonergic_subcluster_summary.csv", index=False)
    enrichment.to_csv(PUBLIC / "serotonergic_marker_enrichment.csv", index=False)
    stability.to_csv(PUBLIC / "serotonergic_resolution_stability.csv", index=False)
    selection_test.to_csv(PUBLIC / "serotonergic_cluster_selection_permutation.csv", index=False)
    libraries.to_csv(PUBLIC / "serotonergic_library_support.csv", index=False)
    markers.to_csv(PUBLIC / "serotonergic_candidate_markers.csv", index=False)
    save_figures(neurons, key, candidate)
    write_report(neurons, candidate, enrichment, stability, libraries, selection_test)

    # Keep the processed object outside public_results because it may be large.
    neurons.write_h5ad(RESULTS / "neurons1_subclustered.h5ad", compression="gzip")
    metadata = {
        "source": str(SOURCE.relative_to(ROOT)),
        "seed": SEED,
        "parent_cluster": TARGET_PARENT_CLUSTER,
        "resolutions": list(RESOLUTIONS),
        "selected_resolution": SELECTED_RESOLUTION,
        "candidate_cluster": candidate,
        "markers": MARKERS,
    }
    (RESULTS / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
