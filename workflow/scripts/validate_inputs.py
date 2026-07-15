"""Fail-fast validation for expression/proteome identifier compatibility."""
import argparse
from pathlib import Path
import anndata as ad
from workflow.scripts.common import fasta_ids, load_config, sha256, write_json

def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); parser.add_argument("--output", required=True); args = parser.parse_args()
    cfg = load_config(args.config); report = {"status": "pass", "species": {}, "errors": []}
    for sid, spec in cfg["species"].items():
        h5ad, fasta = Path(spec["h5ad"]), Path(spec["protein_fasta"])
        entry = {"name": spec["name"], "h5ad": str(h5ad), "protein_fasta": str(fasta)}
        missing = [str(p) for p in (h5ad, fasta) if not p.is_file()]
        if missing:
            report["errors"].append(f"{sid}: missing files: {', '.join(missing)}"); entry["status"] = "missing"; report["species"][sid] = entry; continue
        data = ad.read_h5ad(h5ad, backed="r"); required = [spec["cell_type_key"], spec["sample_key"]]
        absent = [key for key in required if key not in data.obs.columns]
        genes, proteins = set(map(str, data.var_names)), fasta_ids(fasta); overlap = genes & proteins
        entry.update({"cells": int(data.n_obs), "genes": int(data.n_vars), "protein_ids": len(proteins), "matching_ids": len(overlap), "expression_match_fraction": len(overlap) / max(len(genes), 1), "h5ad_sha256": sha256(h5ad), "protein_fasta_sha256": sha256(fasta)})
        if not data.var_names.is_unique: report["errors"].append(f"{sid}: expression gene identifiers are not unique")
        if absent: report["errors"].append(f"{sid}: missing obs columns: {', '.join(absent)}")
        if not overlap: report["errors"].append(f"{sid}: no expression gene IDs match FASTA headers")
        entry["status"] = "pass" if not absent and overlap and data.var_names.is_unique else "fail"; report["species"][sid] = entry; data.file.close()
    maps = Path(cfg["samap"]["maps_directory"])
    if not maps.is_dir() or not any(maps.iterdir()): report["errors"].append(f"homology maps directory missing or empty: {maps}")
    if report["errors"]: report["status"] = "fail"
    write_json(args.output, report)
    if report["errors"]: raise SystemExit("Input validation failed:\n- " + "\n- ".join(report["errors"]))

if __name__ == "__main__": main()
