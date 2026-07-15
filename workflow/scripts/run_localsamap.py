"""Calibrate a saved SAMap graph with blocked permutations and specificity."""
import argparse, itertools, platform
from importlib.metadata import version
from pathlib import Path
import numpy as np
import pandas as pd
from localsamap import benjamini_hochberg, mapping_specificity, permutation_pvalue
from samap.utils import load_samap
from workflow.scripts.common import load_config, write_json

def shuffle_within(values, strata, rng):
    result = np.asarray(values, dtype=object).copy()
    if strata is None: return rng.permutation(result)
    frame = pd.DataFrame({"value": result, "stratum": np.asarray(strata, dtype=object)})
    for indices in frame.groupby("stratum", dropna=False).indices.values(): result[indices] = rng.permutation(result[indices])
    return result

def score_matrix(graph, source_labels, target_labels):
    source_labels, target_labels = source_labels.astype(str), target_labels.astype(str)
    left, right = np.unique(source_labels), np.unique(target_labels); result = np.zeros((left.size, right.size))
    for i, a in enumerate(left):
        rows = source_labels == a
        for j, b in enumerate(right): result[i, j] = float(np.asarray(graph[rows][:, target_labels == b].sum(axis=1)).mean())
    return left, right, result

def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); parser.add_argument("--model", required=True); args = parser.parse_args(); cfg = load_config(args.config)
    params = cfg["localsamap"]; rng = np.random.default_rng(int(cfg["project"]["seed"])); sm = load_samap(args.model); adata = sm.samap.adata
    graph, species, records = adata.obsp["connectivities"].tocsr(), np.asarray(adata.obs["species"]).astype(str), []
    for source, target in itertools.permutations(cfg["species"], 2):
        source_mask, target_mask = species == source, species == target
        source_key, target_key = cfg["species"][source]["cell_type_key"], cfg["species"][target]["cell_type_key"]
        source_labels = np.asarray(sm.sams[source].adata.obs[source_key]).astype(str); target_obs = sm.sams[target].adata.obs; target_labels = np.asarray(target_obs[target_key]).astype(str)
        strata_cols = [x for x in params.get("permutation_strata", []) if x in target_obs]; strata = None if not strata_cols else target_obs[strata_cols].astype(str).agg("|".join, axis=1)
        block = graph[source_mask][:, target_mask]; left, right, observed = score_matrix(block, source_labels, target_labels); null = np.empty((int(params["permutations"]), *observed.shape))
        for r in range(null.shape[0]): _, _, null[r] = score_matrix(block, source_labels, shuffle_within(target_labels, strata, rng))
        for i, a in enumerate(left):
            specificity = mapping_specificity(np.maximum(observed[i], 0) + np.finfo(float).eps)
            for j, b in enumerate(right): records.append({"source_species": source, "target_species": target, "source_cell_type": a, "target_cell_type": b, "score": observed[i, j], "p_value": permutation_pvalue(observed[i, j], null[:, i, j]), "best_target": right[specificity["best_index"]], "specificity_margin": specificity["margin"], "mapping_entropy": specificity["normalized_entropy"], "permutation_strata": ",".join(strata_cols) or "none"})
    results = pd.DataFrame(records); results["q_value"] = benjamini_hochberg(results["p_value"].to_numpy()); results["significant"] = results["q_value"] <= float(params["alpha"])
    out = Path(cfg["project"]["output_root"]) / "localsamap"; out.mkdir(parents=True, exist_ok=True); results.to_csv(out / "mapping_statistics.csv", index=False)
    write_json(out / "run_metadata.json", {"stage": "localsamap_post_analysis", "localsamap_version": version("localsamap"), "python": platform.python_version(), "seed": cfg["project"]["seed"], "parameters": params, "warning": "Validity depends on exchangeability within the recorded permutation strata."})

if __name__ == "__main__": main()
