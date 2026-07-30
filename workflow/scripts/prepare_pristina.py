"""Convert the public Pristina leidyi GSE230505 count matrix to AnnData.

The GEO matrix is transcript by cell; AnnData stores cell by transcript.  This
script deliberately does no quality control or gut-cell selection.  Those are
separate, documented steps so that the downloaded whole-body atlas is never
mistaken for a gut-only dataset.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread


def library_from_barcode(barcode: str) -> str:
    """Return the published library identifier from a GEO barcode when present."""
    parts = barcode.split("_")
    if len(parts) >= 3 and parts[0] == "lib":
        return "_".join(parts[:3])
    return "unknown_library"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--barcodes", required=True, type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    features = pd.read_csv(
        args.features,
        sep="\t",
        header=None,
        names=["transcript_id", "transcript_name", "feature_type"],
        compression="gzip",
        dtype=str,
    )
    barcodes = pd.read_csv(
        args.barcodes,
        sep="\t",
        header=None,
        names=["barcode"],
        compression="gzip",
        dtype=str,
    )
    counts = mmread(args.matrix).tocsr().transpose().tocsr()

    expected = (len(barcodes), len(features))
    if counts.shape != expected:
        raise ValueError(
            f"Matrix shape after transpose is {counts.shape}; expected {expected} "
            "from the barcode and feature tables."
        )
    if features["transcript_id"].duplicated().any():
        raise ValueError("Transcript identifiers are not unique; refusing to alter GEO IDs.")
    if barcodes["barcode"].duplicated().any():
        raise ValueError("Cell barcodes are not unique.")

    obs = pd.DataFrame(index=pd.Index(barcodes["barcode"], name="barcode"))
    obs["library_id"] = obs.index.map(library_from_barcode).astype("category")
    obs["sample_id"] = obs["library_id"]
    obs["species"] = "Pristina leidyi"
    obs["dataset_accession"] = "GSE230505"
    obs["source"] = "whole_body_scRNA_seq"
    # This placeholder is structural only.  It must be replaced with an
    # author-derived or reproducibly calculated annotation before mapping scores
    # are interpreted at cell-type level.
    obs["cell_type"] = "unannotated"

    var = features.set_index("transcript_id", drop=True)
    var.index.name = "transcript_id"
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    adata.uns["provenance"] = {
        "geo_accession": "GSE230505",
        "matrix_orientation_from_geo": "transcript_by_cell",
        "matrix_orientation_in_anndata": "cell_by_transcript",
        "note": "Whole-body atlas imported without QC or gut-cell selection.",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(adata)
    print(f"Wrote {args.output}")
    print(f"Libraries: {obs['library_id'].nunique()}")


if __name__ == "__main__":
    main()
