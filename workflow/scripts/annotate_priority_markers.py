"""Annotate top Pristina marker homologs for focused biological validation.

The script keeps only the top BLAST hit for the top 15 markers from the three
highest-priority SAMap correspondence candidates. It looks up human and
zebrafish Ensembl gene IDs through Ensembl's official REST API, then writes a
small review table. A BLAST hit is a candidate homolog, not a confirmed
ortholog or functional equivalence claim.
"""

from __future__ import annotations

import json
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results_pristina_human_zebrafish/original_samap/pristina_marker_homologs_top3.csv"
OUT = ROOT / "public_results/priority_marker_validation.csv"
PRIORITY_TYPES = ["stomach 2", "anterior/mid-intestine 2", "anterior intestine"]
API_URL = "https://rest.ensembl.org/lookup/id"


def ensembl_lookup(ids: list[str]) -> dict[str, dict]:
    """Return Ensembl's display names and descriptions for up to 1,000 IDs."""
    if not ids:
        return {}
    payload = json.dumps({"ids": ids}).encode("utf-8")
    request = Request(
        API_URL,
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urlopen(request, timeout=60) as response:
        return json.load(response)


def annotation_columns(ids: pd.Series, prefix: str) -> pd.DataFrame:
    lookup = ensembl_lookup(sorted(ids.dropna().astype(str).unique()))
    return pd.DataFrame(
        {
            f"{prefix}_gene_symbol": ids.map(
                lambda gene_id: lookup.get(str(gene_id), {}).get("display_name")
                if pd.notna(gene_id)
                else None
            ),
            f"{prefix}_description": ids.map(
                lambda gene_id: lookup.get(str(gene_id), {}).get("description")
                if pd.notna(gene_id)
                else None
            ),
        }
    )


def main() -> None:
    markers = pd.read_csv(SOURCE, low_memory=False)
    chosen = markers.loc[
        markers["pristina_cell_type"].isin(PRIORITY_TYPES)
        & markers["marker_rank"].le(15)
        & markers["homolog_rank"].eq(1)
    ].copy()
    chosen = chosen.sort_values(["pristina_cell_type", "marker_rank"])
    chosen = chosen.reset_index(drop=True)

    human = annotation_columns(chosen["human_gene"], "human")
    zebrafish = annotation_columns(chosen["zebrafish_gene"], "zebrafish")
    result = pd.concat([chosen, human, zebrafish], axis=1)
    result["annotation_source"] = "Ensembl REST lookup/id"
    result["interpretation_note"] = (
        "Candidate homology from top reciprocal BLAST-table hit; validate before biological label transfer."
    )
    OUT.parent.mkdir(exist_ok=True)
    result.to_csv(OUT, index=False)
    print(f"Wrote: {OUT}")
    print(f"Rows: {len(result)}")
    print(
        result[
            [
                "pristina_cell_type",
                "marker_rank",
                "pristina_gene",
                "human_gene_symbol",
                "zebrafish_gene_symbol",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
