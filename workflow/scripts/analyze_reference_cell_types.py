"""Evaluate professor-suggested Pristina muscle and serotonergic references.

This analysis uses the published Pristina marker/annotation supplements and the
full Pristina atlas. It does not infer cross-species matches when the comparison
atlases lack the required cell types.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "public_results"
REFERENCE = ROOT / "data" / "reference"
ATLAS = ROOT / "data" / "raw" / "pristina" / "pristina_atlas_final.h5ad"
SUPPLEMENT_1 = REFERENCE / "Supplementary Data 1.xlsx"
SUPPLEMENT_3 = REFERENCE / "Supplementary Data 3.xlsx"

SEROTONIN_GENES = {
    "TPH1-like": "PrileiEVm011352t1",
    "DDC-like": "PrileiEVm010940t1",
    "VMAT-like 1": "PrileiEVm008386t1",
    "VMAT-like 2": "PrileiEVm008853t1",
    "DAT-like monoamine transporter": "PrileiEVm006035t1",
}
NEURON_CLUSTERS = ["3", "40", "41", "42", "45", "49", "59"]
MUSCLE_PATTERN = r"myosin|tropomyosin|troponin|paramyosin"


def dense_vector(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray().ravel()
    return np.asarray(matrix).ravel()


def annotation_lookup() -> pd.DataFrame:
    eggnog = pd.read_excel(SUPPLEMENT_1, sheet_name="eggnog")
    eggnog = eggnog.rename(
        columns={"id": "pristina_gene", "Description": "annotation", "Preferred_Name": "preferred_name"}
    )
    return eggnog[["pristina_gene", "annotation", "preferred_name"]]


def muscle_evidence(annotations: pd.DataFrame) -> pd.DataFrame:
    blocks = []
    for cluster, cell_type in [("5", "muscle 1"), ("6", "muscle 2")]:
        markers = pd.read_excel(SUPPLEMENT_3, sheet_name=f"Cluster {cluster}")
        markers = markers.rename(columns={"names": "pristina_gene", "diamond": "diamond_annotation"})
        markers["marker_rank"] = np.arange(1, len(markers) + 1)
        selected = markers[
            markers["diamond_annotation"].astype(str).str.contains(MUSCLE_PATTERN, case=False, regex=True)
        ].head(12).copy()
        selected.insert(0, "candidate", "striated muscle")
        selected.insert(1, "cluster", cluster)
        selected.insert(2, "cell_type", cell_type)
        blocks.append(selected)
    result = pd.concat(blocks, ignore_index=True)
    return result.merge(annotations, on="pristina_gene", how="left")


def serotonin_evidence(adata, annotations: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    cell_clusters = adata.obs["leiden_1.5"].astype(str).to_numpy()
    cell_types = adata.obs["leiden_1.5_names"].astype(str).to_numpy()
    rows = []
    detected = {}
    for marker, gene in SEROTONIN_GENES.items():
        values = dense_vector(adata[:, gene].X)
        detected[marker] = values > 0
        for cluster in NEURON_CLUSTERS:
            mask = cell_clusters == cluster
            rows.append(
                {
                    "candidate": "serotonergic neuron",
                    "cluster": cluster,
                    "cell_type": np.unique(cell_types[mask])[0],
                    "marker": marker,
                    "pristina_gene": gene,
                    "n_cells": int(mask.sum()),
                    "positive_cells": int(np.sum(values[mask] > 0)),
                    "detection_fraction": float(np.mean(values[mask] > 0)),
                    "mean_stored_expression": float(np.mean(values[mask])),
                }
            )
    evidence = pd.DataFrame(rows).merge(annotations, on="pristina_gene", how="left")
    neuron_mask = np.isin(cell_clusters, NEURON_CLUSTERS)
    tph = detected["TPH1-like"] & neuron_mask
    any_vmat = (detected["VMAT-like 1"] | detected["VMAT-like 2"]) & neuron_mask
    ddc = detected["DDC-like"] & neuron_mask
    monoamine = detected["DAT-like monoamine transporter"] & neuron_mask
    counts = {
        "neuronal_cells": int(neuron_mask.sum()),
        "tph_positive_neuronal_cells": int(tph.sum()),
        "tph_and_vmat_positive_cells": int(np.sum(tph & any_vmat)),
        "tph_and_ddc_positive_cells": int(np.sum(tph & ddc)),
        "tph_and_any_supporting_marker_cells": int(np.sum(tph & (any_vmat | ddc | monoamine))),
    }
    return evidence, counts


def candidate_summary(adata, counts: dict[str, int]) -> pd.DataFrame:
    obs = adata.obs
    cluster_sizes = obs["leiden_1.5"].astype(str).value_counts()
    return pd.DataFrame(
        [
            {
                "candidate": "striated muscle",
                "cluster": "5",
                "cell_type": "muscle 1",
                "n_cells": int(cluster_sizes["5"]),
                "evidence_status": "supported reference candidate",
                "key_evidence": "myosin, striated-muscle light chain, troponin, and tropomyosin are top markers",
                "interpretation": "appropriate positive reference after matching comparison-species striated muscle",
            },
            {
                "candidate": "striated muscle",
                "cluster": "6",
                "cell_type": "muscle 2",
                "n_cells": int(cluster_sizes["6"]),
                "evidence_status": "supported reference candidate",
                "key_evidence": "myosin and striated contractile-complex genes dominate the marker list",
                "interpretation": "appropriate positive reference after matching comparison-species striated muscle",
            },
            {
                "candidate": "serotonergic neuron",
                "cluster": "3",
                "cell_type": "neurons 1",
                "n_cells": int(cluster_sizes["3"]),
                "evidence_status": "provisional subpopulation candidate",
                "key_evidence": (
                    f"all {counts['tph_positive_neuronal_cells']} TPH1-like-positive neuronal cells occur in cluster 3; "
                    f"{counts['tph_and_vmat_positive_cells']} also express a VMAT-like transcript and "
                    f"{counts['tph_and_ddc_positive_cells']} also express DDC-like"
                ),
                "interpretation": "likely a rare serotonergic subset within heterogeneous neurons 1; subclustering is required",
            },
        ]
    )


def write_report(summary: pd.DataFrame, counts: dict[str, int]) -> None:
    text = f"""# Professor-suggested reference cell types

## Outcome

- Pristina clusters 5 (`muscle 1`) and 6 (`muscle 2`) are supported striated-muscle reference candidates. Their published marker lists are dominated by myosin, a striated-muscle regulatory light chain, troponin, tropomyosin, and related contractile genes.
- Cluster 3 (`neurons 1`) contains the leading serotonergic candidate, but the evidence identifies a rare subpopulation rather than the entire cluster. All {counts['tph_positive_neuronal_cells']} neuronal cells with detectable `PrileiEVm011352t1` (annotated as TPH1) are in cluster 3. Of those, {counts['tph_and_vmat_positive_cells']} also have detectable VMAT-like expression and {counts['tph_and_ddc_positive_cells']} have detectable DDC-like expression.

## Important limitation

The current human dataset contains intestinal epithelial populations, and the current zebrafish dataset contains gut epithelial, immune, stromal, and pancreatic populations. Neither comparison dataset includes an annotated striated-muscle or neuronal population. Therefore, the existing gut-only SAMap run cannot provide a valid positive cross-species score for these suggested references.

## Next analysis

1. Subcluster the 5,138 cells in Pristina `neurons 1` and test whether TPH1-, VMAT-, and DDC-positive cells form a reproducible group across experimental libraries.
2. Acquire human and zebrafish reference datasets containing striated muscle and serotonergic neurons, with donor/sample metadata retained.
3. Prepare balanced `.h5ad` inputs and rerun SAMap on the expanded reference panel.
4. Use the muscle pairs as high-confidence positives. Keep the serotonin pair provisional until cluster specificity and marker coexpression are independently supported.
5. Construct negative pairs from biologically distinct types only after defining the positive reference populations.

## Sources

- Pristina atlas paper: https://doi.org/10.1038/s41467-024-47401-6
- GEO expression data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230505
- Supplementary Data 1: transcript annotation table.
- Supplementary Data 3: cluster differential-expression marker tables.
- Supplementary Data 4: top pseudobulk marker tables.

SAMap values are similarity scores, not calibrated probabilities. These candidate labels describe evidence available at this stage and should not be treated as confirmed evolutionary identity.
"""
    (PUBLIC / "reference_cell_type_candidates.md").write_text(text)


def main() -> None:
    required = [ATLAS, SUPPLEMENT_1, SUPPLEMENT_3]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input(s): " + ", ".join(missing))
    PUBLIC.mkdir(parents=True, exist_ok=True)
    annotations = annotation_lookup()
    adata = sc.read_h5ad(ATLAS)
    missing_genes = sorted(set(SEROTONIN_GENES.values()) - set(adata.var_names))
    if missing_genes:
        raise ValueError(f"Serotonin candidate genes absent from atlas: {missing_genes}")

    muscle = muscle_evidence(annotations)
    serotonin, counts = serotonin_evidence(adata, annotations)
    summary = candidate_summary(adata, counts)

    muscle.to_csv(PUBLIC / "reference_muscle_marker_evidence.csv", index=False)
    serotonin.to_csv(PUBLIC / "reference_serotonin_marker_evidence.csv", index=False)
    summary.to_csv(PUBLIC / "reference_cell_type_candidates.csv", index=False)
    write_report(summary, counts)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
