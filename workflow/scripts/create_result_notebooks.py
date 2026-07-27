"""Create reader-facing notebooks for compact public SAMap result tables."""

from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS = ROOT / "notebooks"


def markdown(text: str):
    return nbf.v4.new_markdown_cell(dedent(text).strip())


def code(text: str):
    return nbf.v4.new_code_cell(dedent(text).strip())


SETUP = """
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path.cwd().resolve()
if not (ROOT / "public_results").exists():
    ROOT = ROOT.parent
PUBLIC = ROOT / "public_results"

plt.rcParams.update({"figure.dpi": 120, "axes.spines.top": False, "axes.spines.right": False})
"""


def write_notebook(name: str, cells: list) -> None:
    notebook = nbf.v4.new_notebook(cells=cells)
    notebook.metadata = {
        "kernelspec": {"display_name": "samap-gut", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"},
    }
    nbf.write(notebook, NOTEBOOKS / name)


def main() -> None:
    NOTEBOOKS.mkdir(exist_ok=True)

    write_notebook("01_samap_mapping_summary.ipynb", [
        markdown("""
        # SAMap mapping summary

        ## tl;dr

        The official three-species SAMap run identifies several interpretable
        cross-species candidates. The strongest mapping is Pristina
        **anterior/mid-intestine 2 → zebrafish LRE** (score 0.762). Scores are
        similarities, not probabilities or identity calls.
        """),
        markdown("""
        ## Context & Methods

        This notebook reads the compact public export generated after the
        completed SAMap run. `top_score` is the highest cell-type mapping score
        for a Pristina population and a comparison species. `score_margin` is
        the difference between its best and second-best counterpart and is a
        descriptive measure of mapping specificity.
        """),
        code(SETUP),
        markdown("## Data"),
        code("""
        summary = pd.read_csv(PUBLIC / "pristina_cross_species_summary.csv")
        summary.sort_values(["comparison_species", "top_score"], ascending=[True, False]).head(10)
        """),
        markdown("## Results"),
        code("""
        display_columns = [
            "pristina_cell_type", "comparison_species", "top_match",
            "top_score", "second_match", "score_margin", "specificity_class",
        ]
        summary[display_columns].sort_values("top_score", ascending=False)
        """),
        code("""
        plot_data = summary.sort_values("top_score")
        colors = plot_data["comparison_species"].map({"human": "#4C78A8", "zebrafish": "#F58518"})
        labels = plot_data["pristina_cell_type"] + " → " + plot_data["top_match"] + " (" + plot_data["comparison_species"] + ")"

        fig, ax = plt.subplots(figsize=(9, 6))
        ax.barh(labels, plot_data["top_score"], color=colors)
        ax.set_xlabel("SAMap mapping score")
        ax.set_title("Best SAMap counterpart for each Pristina population")
        ax.set_xlim(0, max(plot_data["top_score"]) * 1.12)
        fig.tight_layout()
        plt.show()
        """),
        code("""
        margin_data = summary.sort_values("score_margin")
        colors = margin_data["comparison_species"].map({"human": "#4C78A8", "zebrafish": "#F58518"})
        labels = margin_data["pristina_cell_type"] + " (" + margin_data["comparison_species"] + ")"

        fig, ax = plt.subplots(figsize=(9, 6))
        ax.barh(labels, margin_data["score_margin"], color=colors)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Best-minus-second mapping score")
        ax.set_title("Descriptive specificity of SAMap mapping candidates")
        fig.tight_layout()
        plt.show()
        """),
        markdown("""
        ## Takeaways

        - **Anterior/mid-intestine 2 → zebrafish LRE** has the highest mapping
          score and the largest margin in this comparison.
        - **Stomach 2 → human colonocyte** and **anterior intestine → human
          enterocyte** are useful hypotheses, but marker-level evidence is
          needed before assigning biological labels.
        - Low score margins indicate ambiguity; they should not be presented as
          one-to-one homologous cell types.
        """),
    ])

    write_notebook("02_cross_species_module_validation.ipynb", [
        markdown("""
        # Cross-species module validation

        ## tl;dr

        Translation of reference marker modules through the BLAST maps supports
        all four pre-specified correspondence tests. The zebrafish LRE module
        in Pristina anterior/mid-intestine 2 is the strongest signal
        (permutation z = 47.3; FDR = 0.0005).
        """),
        markdown("""
        ## Context & Methods

        Reference matrices were normalized to 10,000 counts per cell and
        log1p-transformed. For each pre-specified reference type, the analysis
        selected up-regulated markers, translated each marker to its best
        Pristina BLAST hit, and averaged gene-wise z-scored Pristina expression.

        The one-sided permutation p-value is
        `p = (1 + #{T_perm ≥ T_obs}) / (R + 1)`, with `R = 2,000` random
        label sets of the same size as the expected Pristina population.
        """),
        code(SETUP),
        markdown("## Data"),
        code("""
        enrichment = pd.read_csv(PUBLIC / "cross_species_module_enrichment.csv")
        modules = pd.read_csv(PUBLIC / "translated_reference_marker_modules.csv")
        enrichment
        """),
        markdown("## Results"),
        code("""
        labels = (
            enrichment["reference_cell_type"] + " (" +
            enrichment["reference_species"].map({"hs": "human", "dr": "zebrafish"}) + ") → " +
            enrichment["expected_pristina_cell_type"]
        )
        colors = enrichment["reference_species"].map({"hs": "#4C78A8", "dr": "#F58518"})
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.barh(labels, enrichment["permutation_z_score"], color=colors)
        ax.invert_yaxis()
        ax.set_xlabel("Cell-level permutation z-score")
        ax.set_title("Enrichment of translated reference modules in Pristina")
        for i, value in enumerate(enrichment["permutation_z_score"]):
            ax.text(value + 0.6, i, f"{value:.1f}", va="center")
        fig.tight_layout()
        plt.show()
        """),
        code("""
        enrichment[[
            "reference_cell_type", "expected_pristina_cell_type",
            "translated_module_genes", "expected_pristina_cells",
            "permutation_z_score", "one_sided_permutation_p", "bh_fdr",
        ]]
        """),
        code("""
        module_counts = modules.groupby(["test_id", "reference_cell_type", "expected_pristina_cell_type"]).size().reset_index(name="translated_genes")
        module_counts
        """),
        markdown("""
        ## Takeaways

        - The LRE signal is substantially stronger than the other planned tests.
        - The human Paneth, enterocyte, and colonocyte modules provide
          convergent program-level support for their respective candidates.
        - These are **cell-level exploratory** tests. They do not replace
          animal/donor-aware pseudobulk testing, and they do not prove that
          Pristina cells have identical mammalian or zebrafish identities.
        """),
    ])

    write_notebook("03_priority_marker_review.ipynb", [
        markdown("""
        # Priority marker and homolog review

        ## tl;dr

        This notebook reviews the top Pristina markers supporting the three
        priority candidate correspondences. It is intended for choosing a short
        validation panel, not for automatically transferring cell-type names.
        """),
        markdown("""
        ## Context & Methods

        The input table contains differentially expressed Pristina markers and
        their top BLAST-supported human and zebrafish candidate homologs.
        Functional descriptions are Ensembl REST annotations. A BLAST hit is a
        homology candidate, not proof of orthology or conserved function.
        """),
        code(SETUP),
        markdown("## Data"),
        code("""
        markers = pd.read_csv(PUBLIC / "priority_marker_validation.csv")
        markers[[
            "pristina_cell_type", "marker_rank", "pristina_gene",
            "human_gene_symbol", "zebrafish_gene_symbol", "logfoldchanges", "pvals_adj",
        ]].head(12)
        """),
        markdown("## Results"),
        code("""
        top = markers.loc[markers["marker_rank"] <= 10].copy()
        top["label"] = top["pristina_cell_type"] + ": " + top["pristina_gene"]
        top = top.sort_values("logfoldchanges")

        fig, ax = plt.subplots(figsize=(9, 7))
        ax.barh(top["label"], top["logfoldchanges"], color="#4C78A8")
        ax.set_xlabel("Pristina marker log fold-change")
        ax.set_title("Top priority markers by Pristina cell population")
        fig.tight_layout()
        plt.show()
        """),
        code("""
        review = markers[[
            "pristina_cell_type", "marker_rank", "pristina_gene",
            "human_gene_symbol", "human_description",
            "zebrafish_gene_symbol", "zebrafish_description", "interpretation_note",
        ]].sort_values(["pristina_cell_type", "marker_rank"])
        review.head(20)
        """),
        markdown("""
        ## Takeaways

        - The anterior/mid-intestine 2 marker set contains a coherent
          lysosomal/catabolic program, which is consistent with its zebrafish
          LRE mapping hypothesis.
        - The stomach 2 and anterior-intestine candidates require more specific
          marker confirmation because some leading genes are broadly expressed
          metabolic or housekeeping genes.
        - Use the reviewed table to select targets for spatial or in situ
          experiments, while checking each gene's expression pattern directly
          in the Pristina atlas.
        """),
    ])

    write_notebook("04_lre_validation_panel.ipynb", [
        markdown("""
        # LRE-like spatial validation panel

        ## tl;dr

        The primary positive panel for the Pristina anterior/mid-intestine 2
        LRE-like hypothesis contains five translated lysosomal/protein-digestion
        genes: `fuca1.1`, `hexb`, `ctsl.1`, `ctsbb`, and `naga`. The canonical
        zebrafish LRE machinery candidates `cubn`, `dab2`, and `amn` are not
        enriched in this Pristina population, so they are context checks rather
        than primary positive probes.
        """),
        markdown("""
        ## Context & Methods

        The panel begins with zebrafish LRE markers translated to Pristina by
        the precomputed BLAST table. Pristina expression was normalized to
        10,000 counts per cell and log1p transformed. Primary panel candidates
        require high alignment support, positive target-versus-other expression,
        and a detection-fraction difference above 0.05.

        **Key assumption:** BLAST alignment supports a candidate homology
        relationship; it does not prove one-to-one orthology or conserved
        function. The literature rationale for the zebrafish LRE program is
        Park et al. (2019), *Developmental Cell*,
        doi:10.1016/j.devcel.2019.08.001.
        """),
        code(SETUP),
        markdown("## Data"),
        code("""
        panel = pd.read_csv(PUBLIC / "lre_validation_panel.csv")
        context = pd.read_csv(PUBLIC / "lre_mechanistic_context_checks.csv")
        panel[[
            "panel_rank", "gene_symbol", "pristina_gene", "identity",
            "alignment_length", "pristina_log2_fc_target_vs_other",
            "pristina_target_detection_fraction", "pristina_detection_difference",
        ]]
        """),
        markdown("## Results"),
        code("""
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
        plot_data = panel.sort_values("pristina_log2_fc_target_vs_other")
        axes[0].barh(plot_data["gene_symbol"], plot_data["pristina_log2_fc_target_vs_other"], color="#4C78A8")
        axes[0].set_xlabel("Pristina log2 FC: target vs other cells")
        axes[0].set_title("Expression enrichment")
        axes[1].barh(plot_data["gene_symbol"], plot_data["pristina_detection_difference"], color="#F58518")
        axes[1].set_xlabel("Detection-fraction difference")
        axes[1].set_title("Detection specificity")
        fig.suptitle("Primary Pristina LRE-like spatial validation panel", y=1.02)
        fig.tight_layout()
        plt.show()
        """),
        code("""
        context[[
            "gene_symbol", "pristina_gene", "homology_quality",
            "pristina_log2_fc_target_vs_other", "pristina_target_detection_fraction",
            "context_check",
        ]]
        """),
        markdown("""
        ## Takeaways

        - Use several primary panel genes together; a multiplex spatial pattern
          is more informative than any individual gene.
        - The panel tests whether anterior/mid-intestine 2 has a coherent
          lysosomal/protein-digestion program.
        - The absence of `cubn`, `dab2`, or `amn` enrichment means the current
          data do not support simply calling these Pristina cells canonical
          zebrafish LREs. Their function and orthology require independent
          investigation.
        """),
    ])
    print(f"Wrote notebooks to {NOTEBOOKS}")


if __name__ == "__main__":
    main()
