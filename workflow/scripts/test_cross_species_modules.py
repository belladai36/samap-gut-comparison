"""Test reference cell-type marker modules in Pristina gut populations.

This is an exploratory cell-level permutation analysis.  It derives reference
markers from the input AnnData objects, translates those genes with the
precomputed BLAST tables, scores the translated module in Pristina cells, and
tests whether its pre-specified Pristina population has an unusually high
mean score under random relabelling of Pristina cells.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "public_results"
PROCESSED = ROOT / "data" / "processed"
MAPS = ROOT / "data" / "homology" / "maps"
SEED = 20260715
N_PERMUTATIONS = 2_000
TOP_REFERENCE_MARKERS = 50

# Candidate correspondences were pre-specified from the SAMap output.
TESTS = [
    ("human_colonocyte", "hs", "Colonocyte", "stomach 2"),
    ("human_paneth", "hs", "Paneth", "anterior/mid-intestine 2"),
    ("human_enterocyte", "hs", "Enterocyte", "anterior intestine"),
    ("zebrafish_lre", "dr", "LRE", "anterior/mid-intestine 2"),
]


def map_path(species: str) -> Path:
    return MAPS / ("hspr" if species == "hs" else "drpr") / f"{species}_to_pr.txt"


def top_markers(reference: sc.AnnData, cell_type: str) -> pd.DataFrame:
    """Return up-regulated reference markers for one annotated cell type."""
    reference = reference.copy()
    sc.tl.rank_genes_groups(
        reference,
        groupby="cell_type",
        groups=[cell_type],
        reference="rest",
        method="t-test_overestim_var",
        n_genes=TOP_REFERENCE_MARKERS,
        pts=False,
    )
    table = sc.get.rank_genes_groups_df(reference, group=cell_type)
    return table.loc[
        (table["pvals_adj"] < 0.05) & (table["logfoldchanges"] > 0),
        ["names", "logfoldchanges", "pvals_adj"],
    ].rename(columns={"names": "reference_gene"})


def translate_markers(markers: pd.DataFrame, species: str) -> pd.DataFrame:
    """Map each source marker to its best Pristina BLAST target by E-value."""
    columns = [
        "reference_gene", "pristina_gene", "identity", "alignment_length",
        "mismatches", "gap_opens", "query_start", "query_end", "target_start",
        "target_end", "evalue", "bit_score",
    ]
    blast = pd.read_csv(map_path(species), sep="\t", header=None, names=columns)
    blast = blast.sort_values(["reference_gene", "evalue", "bit_score"], ascending=[True, True, False])
    blast = blast.drop_duplicates("reference_gene", keep="first")
    translated = markers.merge(blast, on="reference_gene", how="inner")
    return translated.drop_duplicates("pristina_gene", keep="first")


def z_module_score(matrix, gene_indices: np.ndarray) -> np.ndarray:
    """Mean per-gene z-scored log-expression for a translated gene module."""
    x = matrix[:, gene_indices]
    if sparse.issparse(x):
        means = np.asarray(x.mean(axis=0)).ravel()
        second_moments = np.asarray(x.power(2).mean(axis=0)).ravel()
        variances = np.maximum(second_moments - means**2, 1e-12)
        scores = (x.sum(axis=1).A1 - len(gene_indices) * means.mean())
        # Exact standardized sum, allowing a different variance for each gene.
        scaled = x.multiply(1 / np.sqrt(variances)).sum(axis=1).A1
        return (scaled - np.sum(means / np.sqrt(variances))) / len(gene_indices)
    means = x.mean(axis=0)
    stds = np.maximum(x.std(axis=0), 1e-6)
    return ((x - means) / stds).mean(axis=1)


def permutation_test(scores: np.ndarray, membership: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    """One-sided random-labelling test with the observed label-set size fixed."""
    n_selected = int(membership.sum())
    observed = float(scores[membership].mean())
    null = np.empty(N_PERMUTATIONS, dtype=float)
    for i in range(N_PERMUTATIONS):
        null[i] = rng.choice(scores, size=n_selected, replace=False).mean()
    p_value = (1 + np.count_nonzero(null >= observed)) / (N_PERMUTATIONS + 1)
    z_score = (observed - null.mean()) / null.std(ddof=1)
    return observed, float(z_score), float(p_value)


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    return pd.Series(adjusted, index=p_values.index)


def normalized_log1p(adata: sc.AnnData) -> sc.AnnData:
    """Return a copy normalized to 10,000 counts per cell and log1p transformed."""
    result = adata.copy()
    sc.pp.normalize_total(result, target_sum=10_000)
    sc.pp.log1p(result)
    return result


def main() -> None:
    PUBLIC.mkdir(exist_ok=True)
    pristina = normalized_log1p(sc.read_h5ad(PROCESSED / "pristina_gut_samap.h5ad"))
    references = {
        "hs": normalized_log1p(sc.read_h5ad(PROCESSED / "human_gut_balanced_1000.h5ad")),
        "dr": normalized_log1p(sc.read_h5ad(PROCESSED / "zebrafish_gut.h5ad")),
    }
    rng = np.random.default_rng(SEED)
    rows, translated_rows = [], []

    for test_id, species, source_type, expected_type in TESTS:
        markers = top_markers(references[species], source_type)
        translated = translate_markers(markers, species)
        translated["test_id"] = test_id
        translated["reference_species"] = species
        translated["reference_cell_type"] = source_type
        translated["expected_pristina_cell_type"] = expected_type
        translated_rows.append(translated)

        gene_indices = pristina.var_names.get_indexer(translated["pristina_gene"])
        gene_indices = gene_indices[gene_indices >= 0]
        scores = z_module_score(pristina.X, gene_indices)
        membership = pristina.obs["cell_type"].astype(str).to_numpy() == expected_type
        observed, z_score, p_value = permutation_test(scores, membership, rng)
        rows.append({
            "test_id": test_id,
            "reference_species": species,
            "reference_cell_type": source_type,
            "expected_pristina_cell_type": expected_type,
            "reference_markers_tested": len(markers),
            "translated_module_genes": len(gene_indices),
            "expected_pristina_cells": int(membership.sum()),
            "observed_mean_z_module_score": observed,
            "permutation_z_score": z_score,
            "one_sided_permutation_p": p_value,
            "permutations": N_PERMUTATIONS,
        })

    results = pd.DataFrame(rows)
    results["bh_fdr"] = benjamini_hochberg(results["one_sided_permutation_p"])
    results.to_csv(PUBLIC / "cross_species_module_enrichment.csv", index=False)
    pd.concat(translated_rows, ignore_index=True).to_csv(PUBLIC / "translated_reference_marker_modules.csv", index=False)

    plot_labels = (
        results["reference_cell_type"] + " (" + results["reference_species"].map({"hs": "human", "dr": "zebrafish"})
        + ") → " + results["expected_pristina_cell_type"]
    )
    figure, axis = plt.subplots(figsize=(8.2, 3.8))
    colors = ["#4C78A8" if species == "hs" else "#F58518" for species in results["reference_species"]]
    axis.barh(plot_labels, results["permutation_z_score"], color=colors)
    axis.invert_yaxis()
    axis.set_xlabel("Cell-level permutation z-score")
    axis.set_title("Translated reference modules are enriched in pre-specified Pristina populations")
    axis.text(0.99, -0.25, "2,000 label permutations; exploratory cell-level test", transform=axis.transAxes,
              ha="right", va="top", fontsize=8)
    figure.tight_layout()
    figure.savefig(PUBLIC / "cross_species_module_enrichment.png", dpi=200, bbox_inches="tight")
    plt.close(figure)

    methods = f"""# Cross-species module-enrichment test

This exploratory analysis tests four candidate correspondences selected from the
SAMap mapping results. For each human or zebrafish reference cell type, we
normalized every count matrix to 10,000 counts per cell and log1p-transformed
it. We then derived up-regulated markers with Scanpy's `t-test_overestim_var` comparison
against all other reference cell types (adjusted p < 0.05; positive log fold
change; at most {TOP_REFERENCE_MARKERS} genes). Each source gene was translated
to the best Pristina BLAST hit by minimum E-value. The Pristina score is the
mean gene-wise z-scored expression over the translated module.

For a pre-specified Pristina population of size n, the statistic is its mean
module score. The one-sided p-value uses {N_PERMUTATIONS:,} random label-set
permutations of size n:

`p = (1 + #{{T_perm >= T_obs}}) / ({N_PERMUTATIONS} + 1)`.

BH FDR is applied across the four planned tests. This is a *cell-level,
exploratory* test: its label-exchangeability null does not make cells biological
replicates, and it does not control for donor, batch, or evolutionary dependence.
In particular, the zebrafish dataset has one sample. Candidate identities still
require independent biological validation, preferably donor/animal-aware
pseudobulk and spatial or in situ evidence.
"""
    (PUBLIC / "cross_species_module_enrichment_methods.md").write_text(methods)
    (PUBLIC / "cross_species_module_enrichment_metadata.json").write_text(json.dumps({
        "seed": SEED, "permutations": N_PERMUTATIONS,
        "top_reference_markers": TOP_REFERENCE_MARKERS,
        "tests": [list(test) for test in TESTS],
    }, indent=2) + "\n")
    print(results.to_string(index=False))
    print(f"Wrote: {PUBLIC / 'cross_species_module_enrichment.csv'}")


if __name__ == "__main__":
    main()
