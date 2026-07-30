# Pristina gut cross-species SAMap project

This project compares Pristina gut cell types with human and zebrafish
intestinal cell atlases using SAMap.

The current study status, results, limitations, and next steps are maintained
in the [living project report](PROGRESS.md). Reproducible result notebooks are
in [`notebooks/`](notebooks/).

## Environment

Create the project-local environment (the included `.condarc` keeps both
the environment and package cache inside this project):

```bash
CONDARC="$PWD/.condarc" conda env create \
  --prefix "$PWD/.conda/envs/samap-gut" \
  -f environment.yml
```

Activate it:

```bash
conda activate "$PWD/.conda/envs/samap-gut"
```

Verify the installation:

```bash
python workflow/scripts/check_environment.py
```

Without activating it, run the same check with:

```bash
MPLCONFIGDIR="$PWD/.cache/matplotlib" conda run --no-capture-output \
  --prefix "$PWD/.conda/envs/samap-gut" \
  python -u workflow/scripts/check_environment.py
```

## Comparative workflow

Study-specific inputs remain local and are excluded from version control. Start
from the portable template, then configure paths and metadata keys locally:

```bash
cp config/config.example.yaml config/config.yaml
```

Then run:

```bash
XDG_CACHE_HOME="$PWD/.cache" snakemake -n
XDG_CACHE_HOME="$PWD/.cache" snakemake --cores 4
```

The workflow validates expression/protein identifiers and metadata, then runs
the official SAMap implementation and writes results to
`results/original_samap/`. This public repository does not include an
experimental local-pattern extension; that method should be released and
installed separately before it is used in a reproducible workflow.

## Repository checks

GitHub Actions compiles the public Python scripts on every push and pull
request. It does not run SAMap because the full analysis requires large,
externally obtained expression and sequence files that are intentionally not
distributed here.

## Citation, datasets, and permissions

This repository is proprietary: written permission is required before use,
copying, modification, or redistribution. See [LICENSE](LICENSE).

For the recommended citation for this repository, use [CITATION.cff](CITATION.cff).
For complete attribution of the Pristina, human, and zebrafish datasets, the
SAMap repository and paper, BLAST+, Scanpy, AnnData, and Ensembl resources, see
[REFERENCES.md](REFERENCES.md).
