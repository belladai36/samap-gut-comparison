"""Create the SAMap input containing published Pristina gut cells only."""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    source = sc.read_h5ad(args.input, backed="r")
    required = {"broad_names", "leiden_1.5_names"}
    missing = required.difference(source.obs.columns)
    if missing:
        source.file.close()
        raise KeyError(f"Required annotation columns are missing: {sorted(missing)}")

    mask = source.obs["broad_names"].astype(str).eq("gut").to_numpy()
    gut = source[mask].to_memory()
    source.file.close()

    gut.obs["cell_type"] = gut.obs["leiden_1.5_names"].astype(str).astype("category")
    gut.obs["sample_id"] = (
        gut.obs["Library"].astype(str).astype("category")
        if "Library" in gut.obs
        else "unknown_library"
    )
    gut.obs["species"] = "Pristina leidyi"
    gut.obs["dataset_accession"] = "GSE230505"
    gut.obs["selection_rule"] = "broad_names == 'gut'"

    if "counts" not in gut.layers:
        gut.layers["counts"] = gut.X.copy()
    gut.uns["samap_input"] = {
        "source_atlas": str(args.input),
        "selection_rule": "broad_names == 'gut'",
        "cell_type_key": "leiden_1.5_names copied to obs['cell_type']",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    gut.write_h5ad(args.output, compression="gzip")
    print(gut)
    print(f"Wrote {args.output}")
    print("\nGut cell types:")
    print(gut.obs["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()
