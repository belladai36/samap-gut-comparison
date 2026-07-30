"""Prepare and transcriptome-screen candidate regions for LRE-like probes.

This script intentionally produces *regions*, not vendor-specific oligos.  The
final HCR, smFISH, or RNAscope design needs the assay platform's own rules.
"""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "public_results"
REFERENCE = ROOT / "data" / "reference" / "pristina_transcripts.fna"
# Dedicated screening database.  It is rebuilt from the public transcript
# reference when needed, avoiding the older mapping database's mmap issue.
DATABASE = ROOT / "data" / "homology" / "blastdb" / "probe_screen" / "pristina_transcripts"
REGION_LENGTH = 400
# Retain several non-overlapping choices so every target can have a screened
# backup region even when one region overlaps a related transcript.
REGIONS_PER_GENE = 8


def read_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    fragments: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if identifier in wanted:
                    records[identifier] = "".join(fragments).upper()
                identifier = line[1:].split()[0]
                fragments = []
            elif identifier in wanted:
                fragments.append(line)
    if identifier in wanted:
        records[identifier] = "".join(fragments).upper()
    return records


def exact_kmer_screen(regions: list[dict[str, object]], k: int = 20) -> pd.DataFrame:
    """Find non-self transcripts sharing >=3 exact k-mers with each region.

    This conservative fallback is used only if local BLAST fails.  It assesses
    the short perfect-match risk relevant to hybridization probes; it is not a
    replacement for the final assay vendor's off-target algorithm.
    """
    kmer_to_regions: dict[str, set[str]] = {}
    region_to_gene = {str(row["region_id"]): str(row["pristina_gene"]) for row in regions}
    for row in regions:
        region_id, sequence = str(row["region_id"]), str(row["sequence"])
        for i in range(len(sequence) - k + 1):
            kmer_to_regions.setdefault(sequence[i : i + k], set()).add(region_id)

    hits: list[dict[str, object]] = []
    identifier: str | None = None
    fragments: list[str] = []

    def scan_record(subject_id: str | None, sequence: str) -> None:
        if subject_id is None:
            return
        counts: Counter[str] = Counter()
        for i in range(len(sequence) - k + 1):
            for region_id in kmer_to_regions.get(sequence[i : i + k], ()):
                if subject_id != region_to_gene[region_id]:
                    counts[region_id] += 1
        for region_id, count in counts.items():
            if count >= 3:
                hits.append({"region_id": region_id, "subject_id": subject_id, "exact_kmer_matches": count})

    with REFERENCE.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                scan_record(identifier, "".join(fragments).upper())
                identifier, fragments = line[1:].split()[0], []
            else:
                fragments.append(line)
    scan_record(identifier, "".join(fragments).upper())
    return pd.DataFrame(hits, columns=["region_id", "subject_id", "exact_kmer_matches"])


def shannon_entropy(sequence: str) -> float:
    counts = Counter(sequence)
    n = len(sequence)
    return -sum((count / n) * math.log2(count / n) for count in counts.values() if count)


def pick_regions(gene: str, sequence: str) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    for start in range(0, len(sequence) - REGION_LENGTH + 1, 25):
        region = sequence[start : start + REGION_LENGTH]
        gc = (region.count("G") + region.count("C")) / REGION_LENGTH
        entropy = shannon_entropy(region)
        longest_homopolymer = max(len(run) for run in __import__("re").findall(r"A+|C+|G+|T+", region))
        if "N" not in region and 0.35 <= gc <= 0.65 and entropy >= 1.75 and longest_homopolymer <= 8:
            candidates.append({
                "pristina_gene": gene,
                "start_1_based": start + 1,
                "end_1_based": start + REGION_LENGTH,
                "sequence": region,
                "gc_fraction": gc,
                "shannon_entropy": entropy,
                "longest_homopolymer": longest_homopolymer,
            })
    # Favour balanced GC/complexity and then keep non-overlapping regions.
    candidates.sort(key=lambda x: (abs(float(x["gc_fraction"]) - 0.50), -float(x["shannon_entropy"])))
    selected: list[dict[str, object]] = []
    for candidate in candidates:
        start, end = int(candidate["start_1_based"]), int(candidate["end_1_based"])
        if all(end < int(old["start_1_based"]) or start > int(old["end_1_based"]) for old in selected):
            candidate["region_id"] = f"{gene}|{start}-{end}"
            selected.append(candidate)
        if len(selected) == REGIONS_PER_GENE:
            break
    return selected


def main() -> None:
    panel = pd.read_csv(PUBLIC / "lre_validation_panel.csv")
    genes = panel["pristina_gene"].tolist()
    symbols = dict(zip(panel["pristina_gene"], panel["gene_symbol"]))
    sequences = read_fasta(REFERENCE, set(genes))
    absent = sorted(set(genes) - set(sequences))
    if absent:
        raise ValueError(f"Missing target transcript(s): {', '.join(absent)}")

    regions = [region for gene in genes for region in pick_regions(gene, sequences[gene])]
    for region in regions:
        region["zebrafish_symbol"] = symbols[str(region["pristina_gene"])]
        region["transcript_length"] = len(sequences[str(region["pristina_gene"])])
    table = pd.DataFrame(regions)
    fasta = PUBLIC / "lre_probe_candidate_regions.fasta"
    with fasta.open("w") as handle:
        for row in regions:
            handle.write(f">{row['region_id']}\n{row['sequence']}\n")

    # Any strong non-self transcriptome hit is a warning, not an automatic rejection.
    output = PUBLIC / "lre_probe_region_blastn.tsv"
    blastn = shutil.which("blastn")
    screening_method = "exact_20mer_transcriptome_scan"
    exact_details: pd.DataFrame | None = None
    if blastn and (DATABASE.with_suffix(".ndb").exists() or DATABASE.with_suffix(".nhr").exists()):
        process = subprocess.run([
            blastn, "-query", str(fasta), "-db", str(DATABASE), "-task", "blastn",
            "-evalue", "1e-10", "-max_target_seqs", "100", "-outfmt", "6",
            "-out", str(output), "-num_threads", "4",
        ], check=False)
        if process.returncode == 0:
            screening_method = "blastn_transcriptome_screen"
            hits = pd.read_csv(output, sep="\t", header=None, names=[
                "region_id", "subject_id", "identity", "alignment_length", "mismatches",
                "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score",
            ]) if output.exists() and output.stat().st_size else pd.DataFrame(columns=["region_id", "subject_id", "identity", "alignment_length"])
            hits["query_gene"] = hits["region_id"].str.split("|", regex=False).str[0]
            nonself = hits.loc[(hits["subject_id"] != hits["query_gene"]) & (hits["identity"] >= 80) & (hits["alignment_length"] >= 80)]
            warning = nonself.groupby("region_id").size().rename("strong_nonself_hits").reset_index()
        else:
            # Do not retain incomplete BLAST output from a failed database read.
            output.unlink(missing_ok=True)
            exact_details = exact_kmer_screen(regions)
            warning = exact_details.groupby("region_id").size().rename("strong_nonself_hits").reset_index()
    else:
        exact_details = exact_kmer_screen(regions)
        warning = exact_details.groupby("region_id").size().rename("strong_nonself_hits").reset_index()
    table = table.merge(warning, on="region_id", how="left")
    table["strong_nonself_hits"] = table["strong_nonself_hits"].fillna(0).astype(int)
    table["screening_method"] = screening_method
    table["screening_status"] = table["strong_nonself_hits"].map(lambda n: "review_nonself_hits" if n else "no_strong_nonself_hit")
    table = table[[
        "region_id", "zebrafish_symbol", "pristina_gene", "transcript_length",
        "start_1_based", "end_1_based", "gc_fraction", "shannon_entropy",
        "longest_homopolymer", "strong_nonself_hits", "screening_method", "screening_status", "sequence",
    ]]
    table.to_csv(PUBLIC / "lre_probe_candidate_regions.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    if exact_details is not None:
        exact_details.to_csv(PUBLIC / "lre_probe_region_exact_kmer_hits.csv", index=False)

    clean = table.loc[table["strong_nonself_hits"].eq(0)].copy()
    clean["design_rank"] = clean.groupby("pristina_gene").cumcount() + 1
    recommendations = clean.loc[clean["design_rank"].le(2)].copy()
    recommendations["recommended_role"] = recommendations["design_rank"].map({1: "primary region", 2: "backup region"})
    recommendations.to_csv(PUBLIC / "lre_probe_region_recommendations.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    print(table.drop(columns="sequence").to_string(index=False))
    print("\nRecommended regions:")
    print(recommendations.drop(columns="sequence").to_string(index=False))
    detail_name = "lre_probe_region_exact_kmer_hits.csv" if exact_details is not None else output.name
    print(f"\nWrote {fasta.name}, lre_probe_candidate_regions.csv, lre_probe_region_recommendations.csv, and {detail_name}")


if __name__ == "__main__":
    main()
