# LRE-like probe-region preparation

This document records a **platform-neutral** pre-design screen for the five
Pristina LRE-like panel genes. It provides target regions for a probe vendor or
an assay-specific design tool; it is not a final HCR, smFISH, RNAscope, or ISH
oligo order sheet.

## Recommended primary regions

| Gene/program label | Pristina transcript | Primary region (1-based) | Backup region (1-based) |
|---|---|---:|---:|
| fuca1.1 | `PrileiEVm011523t1` | 151–550 | 1976–2375 |
| hexb | `PrileiEVm008813t1` | 101–500 | 901–1300 |
| ctsl.1 | `PrileiEVm014063t1` | 576–975 | 2526–2925 |
| ctsbb | `PrileiEVm014379t1` | 826–1225 | 251–650 |
| naga | `PrileiEVm012889t1` | 826–1225 | 351–750 |

Each region is 400 nt, contains no ambiguous bases, has 0.39–0.49 GC fraction,
and has no long (>8-nt) homopolymer. The full sequences are in
`lre_probe_candidate_regions.fasta`; the exact coordinates and QC values are
in `lre_probe_region_recommendations.csv`.

## Transcriptome screening

The local BLAST nucleotide database returns a memory-mapping error on this
computer, so the final screen used a conservative exact-20-mer scan over all
29,807 Pristina reference transcripts. A non-self transcript sharing at least
three exact 20-mers with a candidate is flagged for review. All primary and
backup regions above have zero flagged non-self transcripts under this rule.
The detailed flagged hits for non-recommended regions are in
`lre_probe_region_exact_kmer_hits.csv`.

This screen reduces obvious short perfect-match risk, but it does not replace
the final assay vendor's specificity algorithm, isoform annotation check, or
wet-lab positive/negative controls.

## Next decision

Choose the assay platform before generating final oligos:

- **HCR RNA-FISH:** use the 400-nt primary regions to design multiple split
  probe pairs per transcript; suitable for multiplex co-localization.
- **smFISH:** tile many 18–22-nt oligos over each primary region and assign
  spectrally separable fluorophores.
- **Chromogenic ISH/RNAscope:** submit the primary-region sequences to the
  provider for their proprietary probe layout and specificity review.

For the first experiment, start with `fuca1.1` and `hexb`, then add `ctsl.1`
as a third marker. Retain the remaining two targets for a second multiplex or
technical backup.
