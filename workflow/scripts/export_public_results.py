"""Export small, shareable summaries and figures from a completed SAMap run.

This script intentionally exports only derived CSV files and static figures.  It
does not copy raw sequencing data, AnnData objects, BLAST tables, or the large
serialized SAMap object.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results_pristina_human_zebrafish" / "original_samap"
PUBLIC = ROOT / "public_results"


def specificity_class(margin: float) -> str:
    """A descriptive, not inferential, specificity label."""
    if margin >= 0.20:
        return "higher specificity"
    if margin >= 0.05:
        return "moderate specificity"
    return "ambiguous"


def configure_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.7, zorder=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.spines["bottom"].set_color("#7a7a7a")


def set_heading(fig: plt.Figure, title: str, subtitle: str) -> None:
    """Place a non-overlapping title and subtitle above the plotting area."""
    fig.suptitle(title, x=0.135, y=0.985, ha="left", fontsize=15, fontweight="bold")
    fig.text(0.135, 0.946, subtitle, ha="left", fontsize=9.5, color="#4f5b66")


def write_readme() -> None:
    (PUBLIC / "README.md").write_text(
        "# Public SAMap result summary\n\n"
        "This folder contains compact, shareable outputs from the completed "
        "three-species SAMap analysis (Pristina, human, zebrafish). It excludes "
        "raw sequencing data, `.h5ad` expression objects, BLAST databases, and "
        "the serialized SAMap object.\n\n"
        "## Files\n\n"
        "- `pristina_cross_species_summary.csv`: best and second-best cell-type "
        "matches, scores, and score margins.\n"
        "- `pristina_marker_homologs_top3.csv`: Pristina gut marker genes with the "
        "top three BLAST-supported human and zebrafish homolog candidates.\n"
        "- `best_match_scores.png`: best match score for each Pristina gut cell "
        "type against each comparison species.\n"
        "- `match_specificity_margins.png`: difference between the best and "
        "second-best match for each cell type and species.\n\n"
        "## Interpretation\n\n"
        "SAMap scores are cross-species similarity scores, not probabilities and "
        "not proof of one-to-one cell-type identity. The `specificity_class` "
        "column is a descriptive ranking heuristic based only on the top-score "
        "margin: higher specificity (>=0.20), moderate specificity (0.05–<0.20), "
        "or ambiguous (<0.05). Candidate matches should be validated with marker "
        "genes and biological context.\n",
        encoding="utf-8",
    )


def main() -> None:
    PUBLIC.mkdir(exist_ok=True)
    confidence = pd.read_csv(RESULTS / "pristina_match_confidence.csv")
    markers = pd.read_csv(RESULTS / "pristina_marker_homologs_top3.csv", low_memory=False)

    summary = confidence.copy()
    summary["specificity_class"] = summary["score_margin"].map(specificity_class)
    summary = summary.sort_values(
        ["comparison_species", "top_score"], ascending=[True, False]
    )
    summary.to_csv(PUBLIC / "pristina_cross_species_summary.csv", index=False)
    markers.to_csv(PUBLIC / "pristina_marker_homologs_top3.csv", index=False)

    order = (
        summary.groupby("pristina_cell_type")["top_score"].max().sort_values().index.tolist()
    )
    species = ["human", "zebrafish"]
    palette = {"human": "#2f6ea9", "zebrafish": "#c46b28"}
    marker_style = {"human": "o", "zebrafish": "s"}
    y = np.arange(len(order))

    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    offsets = {"human": -0.16, "zebrafish": 0.16}
    for comparison in species:
        subset = (
            summary.loc[summary["comparison_species"].eq(comparison)]
            .set_index("pristina_cell_type")
            .reindex(order)
        )
        ax.scatter(
            subset["top_score"],
            y + offsets[comparison],
            s=48,
            marker=marker_style[comparison],
            color=palette[comparison],
            edgecolor="#25313b",
            linewidth=0.5,
            label=comparison.capitalize(),
            zorder=3,
        )
    ax.set_yticks(y, order)
    ax.set_xlim(0, 0.85)
    ax.set_xlabel("Best SAMap mapping score")
    set_heading(
        fig,
        "Best SAMap match score for each Pristina gut cell type",
        "Top cell-type match against human or zebrafish; scores are similarities, not probabilities.",
    )
    configure_axes(ax)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(PUBLIC / "best_match_scores.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    height = 0.30
    for comparison, offset in [("human", -height / 2), ("zebrafish", height / 2)]:
        subset = (
            summary.loc[summary["comparison_species"].eq(comparison)]
            .set_index("pristina_cell_type")
            .reindex(order)
        )
        ax.barh(
            y + offset,
            subset["score_margin"],
            height=height,
            color=palette[comparison],
            edgecolor="#25313b",
            linewidth=0.4,
            label=comparison.capitalize(),
            zorder=2,
        )
    ax.axvline(0.05, color="#707070", linestyle="--", linewidth=0.9, zorder=1)
    ax.axvline(0.20, color="#707070", linestyle="--", linewidth=0.9, zorder=1)
    ax.text(0.052, len(order) - 0.25, "0.05", fontsize=8, color="#555555")
    ax.text(0.202, len(order) - 0.25, "0.20", fontsize=8, color="#555555")
    ax.set_yticks(y, order)
    ax.set_xlabel("Top-score margin (best minus second-best match)")
    set_heading(
        fig,
        "Specificity of the top SAMap match",
        "Larger margins indicate a clearer ranking within that species; dashed lines are descriptive thresholds.",
    )
    configure_axes(ax)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(PUBLIC / "match_specificity_margins.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    write_readme()
    print(f"Wrote public result bundle to: {PUBLIC}")
    for path in sorted(PUBLIC.iterdir()):
        print(f"- {path.name}")


if __name__ == "__main__":
    main()
