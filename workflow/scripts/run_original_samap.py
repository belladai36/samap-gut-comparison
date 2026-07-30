"""Run the official SAMap workflow and write compact result summaries."""

from __future__ import annotations

import argparse
import platform
from importlib.metadata import version
from pathlib import Path

import anndata as ad
from samap import SAMAP, get_mapping_scores
from samap.utils import save_samap

from workflow.scripts.common import load_config, write_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    # AnnData 0.11+ requires an explicit opt-in to serialize nullable strings.
    # SAMap may create these while assembling the combined object.
    ad.settings.allow_write_nullable_strings = True

    cfg = load_config(args.config)
    output = Path(cfg["project"]["output_root"]) / "original_samap"
    output.mkdir(parents=True, exist_ok=True)
    species = cfg["species"]
    params = cfg["samap"]
    inputs = {species_id: spec["h5ad"] for species_id, spec in species.items()}
    keys = {species_id: spec["cell_type_key"] for species_id, spec in species.items()}

    samap = SAMAP(
        sams=inputs,
        f_maps=params["maps_directory"],
        keys=keys,
        eval_thr=float(params["eval_threshold"]),
        save_processed=True,
    )
    samap.run(
        n_iterations=int(params["n_iterations"]),
        cross_species_k=int(params["cross_species_k"]),
        n_gene_chunks=int(params["n_gene_chunks"]),
        ncpus=int(params["ncpus"]),
        umap=bool(params["umap"]),
    )

    save_samap(samap, str(output / "samap.pkl"))
    _, scores = get_mapping_scores(samap, keys=keys)
    scores.to_csv(output / "mapping_scores.csv")
    samap.samap.adata.write_h5ad(output / "combined_samap.h5ad")
    write_json(
        output / "run_metadata.json",
        {
            "stage": "official_samap",
            "samap_version": version("sc-samap"),
            "python": platform.python_version(),
            "parameters": params,
            "species": list(species),
            "cell_type_keys": keys,
        },
    )


if __name__ == "__main__":
    main()
