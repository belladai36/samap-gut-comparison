"""Run official SAMap without LocalSAMap modifications."""
import argparse, platform
from importlib.metadata import version
from pathlib import Path
from samap import SAMAP, get_mapping_scores
from samap.utils import save_samap
from workflow.scripts.common import load_config, write_json

def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); args = parser.parse_args(); cfg = load_config(args.config)
    out = Path(cfg["project"]["output_root"]) / "original_samap"; out.mkdir(parents=True, exist_ok=True)
    species, params = cfg["species"], cfg["samap"]
    inputs = {sid: spec["h5ad"] for sid, spec in species.items()}; keys = {sid: spec["cell_type_key"] for sid, spec in species.items()}
    sm = SAMAP(sams=inputs, f_maps=params["maps_directory"], keys=keys, eval_thr=float(params["eval_threshold"]), save_processed=True)
    sm.run(n_iterations=int(params["n_iterations"]), cross_species_k=int(params["cross_species_k"]), n_gene_chunks=int(params["n_gene_chunks"]), ncpus=int(params["ncpus"]), umap=bool(params["umap"]))
    save_samap(sm, str(out / "samap.pkl")); _, scores = get_mapping_scores(sm, keys=keys); scores.to_csv(out / "mapping_scores.csv"); sm.samap.adata.write_h5ad(out / "combined_samap.h5ad")
    write_json(out / "run_metadata.json", {"stage": "official_samap", "samap_version": version("sc-samap"), "python": platform.python_version(), "parameters": params, "species": list(species), "cell_type_keys": keys})

if __name__ == "__main__": main()
