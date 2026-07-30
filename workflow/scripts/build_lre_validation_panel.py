"""Build a cautious, data-supported Pristina LRE-like validation panel.

The panel is intended for spatial-expression validation. It does not transfer a
zebrafish cell-type label; it prioritizes Pristina genes with both an LRE
program connection and enriched expression in anterior/mid-intestine 2.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "public_results"
PRISTINA = ROOT / "data" / "processed" / "pristina_gut_samap.h5ad"
ZEBRAFISH = ROOT / "data" / "processed" / "zebrafish_gut.h5ad"
TARGET = "anterior/mid-intestine 2"

# Park et al. (2019) identifies Cubn/Amn/Dab2 as LRE endocytic machinery.
CORE_ANCHORS = ["cubn", "dab2", "amn"]
LYSOSOMAL_SUPPORT = ["fuca1.1", "hexb", "ctsl.1", "ctsz", "ctsa", "ctsbb", "lgmn", "naga", "scpep1", "ifi30"]


def column_mean(matrix: sparse.spmatrix | np.ndarray) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.mean(axis=0)).ravel()
    return np.asarray(matrix.mean(axis=0)).ravel()


def fraction_positive(matrix: sparse.spmatrix | np.ndarray) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray((matrix > 0).mean(axis=0)).ravel()
    return np.asarray((matrix > 0).mean(axis=0)).ravel()


def main() -> None:
    modules = pd.read_csv(PUBLIC / "translated_reference_marker_modules.csv")
    lre = modules.loc[modules["test_id"].eq("zebrafish_lre")].copy()

    zebrafish = sc.read_h5ad(ZEBRAFISH, backed="r")
    symbols = zebrafish.var[["gene_symbol"]].reset_index().rename(columns={"ensembl_gene_id": "reference_gene"})
    zebrafish.file.close()
    lre = lre.merge(symbols, on="reference_gene", how="left")

    adata = sc.read_h5ad(PRISTINA)
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    target = adata.obs["cell_type"].astype(str).eq(TARGET).to_numpy()
    indices = adata.var_names.get_indexer(lre["pristina_gene"])
    lre["available_in_pristina_matrix"] = indices >= 0
    lre = lre.loc[lre["available_in_pristina_matrix"]].copy()
    indices = adata.var_names.get_indexer(lre["pristina_gene"])

    x = adata.X[:, indices]
    target_mean = column_mean(x[target])
    other_mean = column_mean(x[~target])
    target_pct = fraction_positive(x[target])
    other_pct = fraction_positive(x[~target])
    lre["pristina_target_mean_log1p"] = target_mean
    lre["pristina_other_mean_log1p"] = other_mean
    lre["pristina_log2_fc_target_vs_other"] = np.log2((target_mean + 0.1) / (other_mean + 0.1))
    lre["pristina_target_detection_fraction"] = target_pct
    lre["pristina_other_detection_fraction"] = other_pct
    lre["pristina_detection_difference"] = target_pct - other_pct
    lre["homology_quality"] = np.where(
        (lre["identity"] >= 35) & (lre["alignment_length"] >= 250), "high", "moderate"
    )

    lre["candidate_class"] = "other LRE-module gene"
    lre.loc[lre["gene_symbol"].isin(CORE_ANCHORS), "candidate_class"] = "core LRE endocytic anchor"
    lre.loc[lre["gene_symbol"].isin(LYSOSOMAL_SUPPORT), "candidate_class"] = "lysosomal/protein-digestion support"

    support = lre.loc[lre["gene_symbol"].isin(LYSOSOMAL_SUPPORT)].copy()
    support = support.loc[
        support["homology_quality"].eq("high")
        & support["pristina_log2_fc_target_vs_other"].gt(0)
        & support["pristina_detection_difference"].gt(0.05)
    ].copy()
    support["support_priority"] = (
        support["pristina_log2_fc_target_vs_other"].rank(pct=True)
        + support["pristina_detection_difference"].rank(pct=True)
        + support["identity"].rank(pct=True)
    )
    support = support.sort_values("support_priority", ascending=False).head(5)
    panel = support.copy()
    panel["panel_rank"] = range(1, len(panel) + 1)
    panel["recommended_use"] = "Primary positive-panel probe: test co-localization with the other panel genes."
    panel["interpretation_caveat"] = (
        "Candidate homolog/program marker only; validate spatial localization and co-expression before assigning LRE-like identity."
    )
    panel.to_csv(PUBLIC / "lre_validation_panel.csv", index=False)

    plot_data = panel.sort_values("pristina_log2_fc_target_vs_other")
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    axes[0].barh(plot_data["gene_symbol"], plot_data["pristina_log2_fc_target_vs_other"], color="#4C78A8")
    axes[0].set_xlabel("Pristina log2 FC: target vs other cells")
    axes[0].set_title("Expression enrichment")
    axes[1].barh(plot_data["gene_symbol"], plot_data["pristina_detection_difference"], color="#F58518")
    axes[1].set_xlabel("Detection-fraction difference")
    axes[1].set_title("Detection specificity")
    figure.suptitle("Primary Pristina LRE-like spatial validation panel", y=1.02)
    figure.tight_layout()
    figure.savefig(PUBLIC / "lre_validation_panel.png", dpi=200, bbox_inches="tight")
    plt.close(figure)

    anchors = lre.loc[lre["gene_symbol"].isin(CORE_ANCHORS)].copy()
    anchors["context_check"] = np.where(
        (anchors["pristina_log2_fc_target_vs_other"] > 0)
        & (anchors["pristina_target_detection_fraction"] >= 0.05),
        "Could be included as a mechanistic anchor after independent orthology review.",
        "Not a primary positive probe: not enriched/detected in this Pristina population.",
    )
    anchors.to_csv(PUBLIC / "lre_mechanistic_context_checks.csv", index=False)

    source = "Park et al. 2019, Developmental Cell, doi:10.1016/j.devcel.2019.08.001"
    lines = [
        "# Pristina LRE-like spatial validation panel",
        "",
        "This panel prioritizes genes for testing the hypothesis that Pristina anterior/mid-intestine 2 contains an LRE-like program.",
        "It is not a definitive cell-type annotation.",
        "",
        "## Selection logic",
        "",
        f"- Target population: `{TARGET}`.",
        "- Primary positive panel: lysosomal/protein-digestion LRE markers with high BLAST alignment support, positive Pristina target-versus-other expression, and >0.05 detection-fraction difference.",
        "- Context checks: zebrafish `cubn`, `dab2`, and `amn` are retained separately because they form LRE-associated receptor-mediated endocytic machinery, but they are not primary positive probes unless they are enriched in Pristina.",
        "- Pristina expression is library-size normalized to 10,000 counts per cell and log1p transformed before scoring.",
        f"- Biological rationale source: {source}.",
        "",
        "## Recommended experiment",
        "",
        "Use a multiplex panel of several primary positive genes. Confirm that signals co-localize in anterior/mid-intestine 2, then compare the localization with established gut-region markers. The canonical zebrafish machinery genes should be interpreted as separate context checks, not assumed to be positive Pristina markers. A negative or discordant result should be interpreted as evidence against, not merely absence of proof for, the LRE-like hypothesis.",
        "",
        "## Panel",
        "",
    ]
    display = panel[[
        "panel_rank", "gene_symbol", "pristina_gene", "candidate_class", "homology_quality",
        "identity", "alignment_length", "pristina_log2_fc_target_vs_other",
        "pristina_target_detection_fraction", "pristina_detection_difference", "recommended_use",
    ]].copy()
    lines += [
        "| Rank | Zebrafish symbol | Pristina gene | Class | Homology quality | Identity (%) | Alignment length | Pristina log2 FC | Target detection | Detection difference |",
        "|---:|---|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in display.itertuples(index=False):
        lines.append(
            f"| {row.panel_rank} | {row.gene_symbol} | {row.pristina_gene} | {row.candidate_class} | "
            f"{row.homology_quality} | {row.identity:.3f} | {row.alignment_length} | "
            f"{row.pristina_log2_fc_target_vs_other:.3f} | {row.pristina_target_detection_fraction:.3f} | "
            f"{row.pristina_detection_difference:.3f} |"
        )
    lines += [
        "",
        "## Caveats",
        "",
        "- BLAST-supported homology is not proof of one-to-one orthology or conserved function.",
        "- The expression statistics are single-cell descriptive quantities, not replicate-level inference.",
        "- `cubn`, `dab2`, and `amn` are reported separately in `lre_mechanistic_context_checks.csv`; they are not used as primary positive probes because their current Pristina expression does not support that use.",
    ]
    (PUBLIC / "lre_validation_panel.md").write_text("\n".join(lines) + "\n")
    print(display.to_string(index=False))
    print("Wrote:", PUBLIC / "lre_validation_panel.csv")


if __name__ == "__main__":
    main()
