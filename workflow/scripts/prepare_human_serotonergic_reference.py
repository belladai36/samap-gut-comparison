"""Convert GSE76381 human ventral-midbrain molecule counts to AnnData.

The GEO CEF file stores four metadata rows followed by a dense gene-by-cell
count matrix.  SAMap requires gene identifiers that agree with the BLAST maps,
so gene symbols are translated to unique Ensembl identifiers using an existing
project human AnnData object as the reference vocabulary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def read_cef(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        cef = next(reader)
        cell_ids = next(reader)
        cell_types = next(reader)
        timepoints = next(reader)

    if cef[:2] != ["CEF", "0"]:
        raise ValueError(f"Unexpected CEF header in {path}")
    if cell_ids[1] != "Cell_ID" or cell_types[1] != "Cell_type":
        raise ValueError("Required Cell_ID/Cell_type metadata rows are missing")
    if timepoints[1] != "Timepoint":
        raise ValueError("Required Timepoint metadata row is missing")

    obs = pd.DataFrame(
        {
            "cell_type": cell_types[2:],
            "timepoint": timepoints[2:],
            "sample_id": [cell_id.split("_")[0] for cell_id in cell_ids[2:]],
        },
        index=pd.Index(cell_ids[2:], name="cell_id"),
    )

    table = pd.read_csv(path, sep="\t", skiprows=4, index_col=0, low_memory=False)
    if table.columns[0].startswith("Unnamed"):
        table = table.iloc[:, 1:]
    if table.shape[1] != len(cell_ids[2:]):
        raise ValueError("Count-matrix width does not match the CEF Cell_ID row")
    # The CEF matrix-header row contains only ``Gene``; recover the cell names
    # from the explicit Cell_ID metadata row above it.
    table.columns = cell_ids[2:]
    if not table.columns.equals(obs.index):
        raise ValueError("Count-matrix columns do not match the CEF Cell_ID row")
    if table.index.has_duplicates:
        raise ValueError("CEF gene symbols are not unique")
    return table, obs


def unique_symbol_to_ensembl(reference_h5ad: Path) -> pd.Series:
    reference = ad.read_h5ad(reference_h5ad, backed="r")
    if "gene_symbol" not in reference.var:
        reference.file.close()
        raise ValueError(f"{reference_h5ad} lacks var['gene_symbol']")
    mapping = pd.DataFrame(
        {
            "symbol": reference.var["gene_symbol"].astype(str).to_numpy(),
            "ensembl": reference.var_names.astype(str),
        }
    )
    reference.file.close()
    mapping = mapping[mapping["symbol"].ne("") & mapping["symbol"].ne("nan")]
    mapping = mapping[~mapping["symbol"].duplicated(keep=False)]
    return mapping.set_index("symbol")["ensembl"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cef", type=Path, required=True)
    parser.add_argument("--gene-reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    counts, obs = read_cef(args.cef)
    symbol_map = unique_symbol_to_ensembl(args.gene_reference)
    matched = counts.index.intersection(symbol_map.index, sort=False)
    if len(matched) == 0:
        raise ValueError("No CEF gene symbols matched the Ensembl reference")

    mapped_counts = counts.loc[matched]
    matrix = sparse.csr_matrix(mapped_counts.to_numpy(dtype=np.float32).T)
    var = pd.DataFrame(
        {"gene_symbol": matched.astype(str)},
        index=pd.Index(symbol_map.loc[matched].to_numpy(), name="ensembl_gene_id"),
    )
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    result.uns["source"] = {
        "accession": "GSE76381",
        "file": args.cef.name,
        "organism": "Homo sapiens",
        "tissue": "fetal ventral midbrain",
        "technology": "STRT-seq",
        "gene_id_translation": "unique exact gene-symbol matches to project human reference",
        "source_gene_count": int(counts.shape[0]),
        "matched_gene_count": int(result.n_vars),
    }
    result.layers["counts"] = result.X.copy()
    result.var_names_make_unique()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(args.output, compression="gzip")

    print(f"Wrote: {args.output}")
    print(f"Shape: {result.n_obs:,} cells x {result.n_vars:,} matched genes")
    print(f"Gene coverage: {result.n_vars / counts.shape[0]:.1%}")
    print("Cell types:")
    print(result.obs["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()
