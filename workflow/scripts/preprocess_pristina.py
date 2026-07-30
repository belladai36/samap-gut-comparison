"""Reproduce the main preprocessing choices in the Pristina atlas notebook.

Input is the unmodified-count AnnData written by ``prepare_pristina.py``.
The workflow keeps raw counts in ``layers['counts']`` and writes a new file,
so it is safe to compare alternative thresholds later.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260716)
    args = parser.parse_args()

    sc.settings.verbosity = 2
    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Values reproduce the public atlas preprocessing notebook.
    sc.pp.filter_cells(adata, min_counts=50)
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, max_counts=1_000_000)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    adata = adata[
        (adata.obs["n_genes_by_counts"] < 700)
        & (adata.obs["total_counts"] < 900)
    ].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=18_000)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=150, random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=45, n_pcs=105, random_state=args.seed)
    sc.tl.umap(adata, min_dist=0.5, spread=1, alpha=1, gamma=1.0, random_state=args.seed)
    for resolution in (0.5, 1.0, 1.5, 2.0):
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=f"leiden_{resolution:g}",
            random_state=args.seed,
        )

    adata.uns["preprocessing"] = {
        "source": "GSE230505 author notebook 3 preprocessing",
        "min_counts": 50,
        "min_genes": 50,
        "max_detected_genes_exclusive": 700,
        "max_total_counts_exclusive": 900,
        "hvg_count": 18_000,
        "neighbors": 45,
        "pcs": 105,
        "seed": args.seed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(adata)
    print(f"Wrote {args.output}")
    for key in ["leiden_0.5", "leiden_1", "leiden_1.5", "leiden_2"]:
        print(f"\n{key}")
        print(adata.obs[key].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
