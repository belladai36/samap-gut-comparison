"""Create the reader-facing serotonergic subclustering analysis notebook."""

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "notebooks/05_serotonergic_neuron_subclustering.ipynb"


def code(source: str):
    return nbf.v4.new_code_cell(source.strip())


def markdown(source: str):
    return nbf.v4.new_markdown_cell(source.strip())


def main() -> None:
    notebook = nbf.v4.new_notebook()
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3 (samap-gut)",
        "language": "python",
        "name": "python3",
    }
    notebook["metadata"]["language_info"] = {"name": "python", "version": "3.12"}
    notebook["cells"] = [
        markdown(
            """
# Pristina serotonergic-neuron subclustering

## tl;dr

Unsupervised reclustering of 5,138 Pristina `neurons 1` cells identifies a stable 226-cell community containing 30 of 46 TPH1-like-positive neurons. TPH1-like detection is 13.3% in this community versus 0.33% elsewhere, and a 10,000-permutation cluster-selection test gives empirical P = 0.0001. The community is present in all nine libraries; TPH1-positive members occur in seven libraries from two experiments. This supports a **provisional serotonergic-like subpopulation**, not a confirmed conserved cell type.
"""
        ),
        markdown(
            """
## Context & Methods

The source is `data/raw/pristina/pristina_atlas_final.h5ad`. Cells annotated as `neurons 1` (parent Leiden cluster 3) were reclustered from raw counts using 3,000 highly variable genes, PCA, a 20-neighbor graph, and Leiden resolutions 0.2–0.8. TPH1/DDC/VMAT genes were **not forced into the clustering feature set** and were evaluated only after clustering.

### Key assumptions

- Detectable TPH1-like expression is treated as evidence for serotonin synthesis potential, not definitive identity.
- Transcript dropout makes marker-negative cells ambiguous.
- Resolution 0.6 was selected as a fixed reporting resolution; the corresponding most enriched cluster was corrected for post-selection using permutations.
- Cross-species identity requires a later SAMap comparison with suitable vertebrate neuronal references.

Reproduce the full calculation with:

```bash
python -m workflow.scripts.subcluster_serotonergic_neurons
```
"""
        ),
        markdown("## Data"),
        code(
            """
from pathlib import Path
import pandas as pd
from IPython.display import Image, display

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
PUBLIC = ROOT / "public_results"

summary = pd.read_csv(PUBLIC / "serotonergic_subcluster_summary.csv")
enrichment = pd.read_csv(PUBLIC / "serotonergic_marker_enrichment.csv")
stability = pd.read_csv(PUBLIC / "serotonergic_resolution_stability.csv")
libraries = pd.read_csv(PUBLIC / "serotonergic_library_support.csv")
selection = pd.read_csv(PUBLIC / "serotonergic_cluster_selection_permutation.csv")
markers = pd.read_csv(PUBLIC / "serotonergic_candidate_markers.csv")

assert len(summary) == 1
assert summary.loc[0, "neurons1_cells"] == 5138
assert selection.loc[0, "n_permutations"] == 10000
summary
"""
        ),
        markdown("## Results\n\n### 1. Marker enrichment"),
        code(
            """
enrichment[[
    "marker", "candidate_positive", "candidate_total",
    "candidate_detection_fraction", "other_detection_fraction",
    "odds_ratio", "fdr_bh"
]].round(5)
"""
        ),
        code(
            """
display(Image(filename=str(PUBLIC / "serotonergic_marker_enrichment.png"), width=850))
"""
        ),
        markdown("### 2. Unsupervised clusters and marker localization"),
        code(
            """
display(Image(filename=str(PUBLIC / "serotonergic_neuron_subclustering_umap.png"), width=1000))
"""
        ),
        markdown("### 3. Resolution stability and selection-corrected inference"),
        code(
            """
display(stability.round(5))
display(selection.round(6))
"""
        ),
        markdown("### 4. Support across experimental libraries"),
        code("libraries"),
        markdown("### 5. Leading differential markers"),
        code(
            """
markers[[
    "names", "logfoldchanges", "pvals_adj", "pct_nz_group",
    "Preferred_name", "Description"
]].head(15)
"""
        ),
        markdown(
            """
## Takeaways

1. The candidate is compact and reproducible at the cell-community level: 194–238 cells across the tested resolutions, with Jaccard similarity of 0.82–0.94 to the resolution-0.6 group.
2. TPH1-like and VMAT-like expression are strongly enriched, and the cluster-search-corrected empirical P value is 0.0001.
3. Only 30/46 TPH1-positive neurons fall inside the candidate, and only 13.3% of candidate cells have detectable TPH1-like transcript. This may reflect dropout, biological heterogeneity, or imperfect clustering.
4. TPH1-positive candidate cells occur in seven libraries but only two of three experiments; the smallest experiment contains just 130 `neurons 1` cells and no detected TPH1-like cells.
5. The evidence is strong enough to carry the group forward as a **provisional serotonergic-like Pristina query**, but cross-species validation and independent spatial or molecular confirmation remain necessary.
"""
        ),
    ]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(notebook, OUTPUT)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
