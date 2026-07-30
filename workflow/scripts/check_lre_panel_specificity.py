"""Confirm cluster-level specificity of the primary Pristina LRE-like panel."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "public_results"
PRISTINA = ROOT / "data" / "processed" / "pristina_gut_samap.h5ad"
TARGET = "anterior/mid-intestine 2"


def mean(matrix):
    return np.asarray(matrix.mean(axis=0)).ravel() if sparse.issparse(matrix) else np.asarray(matrix.mean(axis=0)).ravel()


def fraction_positive(matrix):
    return np.asarray((matrix > 0).mean(axis=0)).ravel() if sparse.issparse(matrix) else np.asarray((matrix > 0).mean(axis=0)).ravel()


def main() -> None:
    panel = pd.read_csv(PUBLIC / "lre_validation_panel.csv")
    genes = panel["pristina_gene"].tolist()
    labels = panel["gene_symbol"].tolist()
    adata = sc.read_h5ad(PRISTINA)
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    cell_types = sorted(adata.obs["cell_type"].astype(str).unique())
    gene_indices = adata.var_names.get_indexer(genes)
    if (gene_indices < 0).any():
        raise ValueError("A primary panel gene is absent from the Pristina expression matrix.")

    rows = []
    for cell_type in cell_types:
        in_group = adata.obs["cell_type"].astype(str).to_numpy() == cell_type
        # scipy sparse matrices treat two index arrays as paired coordinates;
        # np.ix_ requests the intended Cartesian cell-by-gene submatrix.
        x_in = adata.X[np.ix_(np.flatnonzero(in_group), gene_indices)]
        x_out = adata.X[np.ix_(np.flatnonzero(~in_group), gene_indices)]
        means_in, means_out = mean(x_in), mean(x_out)
        detected_in = fraction_positive(x_in)
        detected_out = fraction_positive(x_out)
        for gene, symbol, mu_in, mu_out, pct_in, pct_out in zip(genes, labels, means_in, means_out, detected_in, detected_out):
            rows.append({
                "cell_type": cell_type,
                "pristina_gene": gene,
                "zebrafish_symbol": symbol,
                "n_cells": int(in_group.sum()),
                "mean_log1p_expression": mu_in,
                "detection_fraction": pct_in,
                "log2_fc_vs_other": np.log2((mu_in + 0.1) / (mu_out + 0.1)),
                "detection_difference_vs_other": pct_in - pct_out,
            })
    specificity = pd.DataFrame(rows)
    specificity["mean_expression_rank"] = specificity.groupby("pristina_gene")["mean_log1p_expression"].rank(method="min", ascending=False).astype(int)
    specificity["is_target_population"] = specificity["cell_type"].eq(TARGET)
    specificity.to_csv(PUBLIC / "lre_panel_cluster_specificity.csv", index=False)

    plot_data = specificity.pivot(index="zebrafish_symbol", columns="cell_type", values="mean_log1p_expression").loc[labels, cell_types]
    size_data = specificity.pivot(index="zebrafish_symbol", columns="cell_type", values="detection_fraction").loc[labels, cell_types]
    figure, axis = plt.subplots(figsize=(12, 4.2))
    y, x = np.indices(plot_data.shape)
    scatter = axis.scatter(
        x.ravel(), y.ravel(), s=20 + 480 * size_data.to_numpy().ravel(),
        c=plot_data.to_numpy().ravel(), cmap="Blues", edgecolor="#333333", linewidth=0.25,
    )
    target_col = cell_types.index(TARGET)
    axis.axvspan(target_col - 0.5, target_col + 0.5, color="#F58518", alpha=0.12, zorder=-1)
    axis.set_xticks(range(len(cell_types)), cell_types, rotation=35, ha="right")
    axis.set_yticks(range(len(labels)), labels)
    axis.set_title("Primary LRE-like panel expression across Pristina gut populations")
    axis.set_xlabel("Pristina cell population; orange band = target population")
    axis.set_ylabel("Zebrafish LRE-marker symbol")
    colorbar = figure.colorbar(scatter, ax=axis, pad=0.01)
    colorbar.set_label("Mean normalized log1p expression")
    figure.tight_layout()
    figure.savefig(PUBLIC / "lre_panel_cluster_specificity.png", dpi=200, bbox_inches="tight")
    plt.close(figure)

    figure, axes = plt.subplots(1, len(genes), figsize=(16, 3.4), constrained_layout=True)
    coordinates = adata.obsm["X_umap"]
    values = adata.X[:, gene_indices]
    if sparse.issparse(values):
        values = values.toarray()
    for axis, gene, symbol, value in zip(axes, genes, labels, values.T):
        order = np.argsort(value)
        image = axis.scatter(coordinates[order, 0], coordinates[order, 1], c=value[order], s=2, cmap="Blues", linewidths=0)
        axis.set_title(f"{symbol}\n{gene}", fontsize=9)
        axis.set_xticks([]); axis.set_yticks([])
        figure.colorbar(image, ax=axis, fraction=0.046, pad=0.03)
    figure.suptitle("UMAP expression of primary LRE-like panel genes in Pristina gut", y=1.03)
    figure.savefig(PUBLIC / "lre_panel_umap_expression.png", dpi=200, bbox_inches="tight")
    plt.close(figure)

    target_rows = specificity.loc[specificity["is_target_population"], [
        "zebrafish_symbol", "pristina_gene", "mean_log1p_expression", "detection_fraction",
        "log2_fc_vs_other", "detection_difference_vs_other", "mean_expression_rank",
    ]]
    target_rows.to_csv(PUBLIC / "lre_panel_target_summary.csv", index=False)
    print(target_rows.to_string(index=False))


if __name__ == "__main__":
    main()
