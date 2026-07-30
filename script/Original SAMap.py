"""Validate configured inputs and launch the official three-species SAMap run."""

from __future__ import annotations

import itertools
import subprocess
import sys
from pathlib import Path

from workflow.scripts.common import load_config


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "config.yaml"


def required_inputs(config: dict) -> list[Path]:
    required = [ROOT / spec["h5ad"] for spec in config["species"].values()]
    maps_root = ROOT / config["samap"]["maps_directory"]
    for source, target in itertools.combinations(sorted(config["species"]), 2):
        pair_dirs = [
            path for path in maps_root.iterdir()
            if path.is_dir() and source in path.name and target in path.name
        ] if maps_root.is_dir() else []
        if not pair_dirs:
            required.append(maps_root / f"{source}{target}")
            continue
        pair = pair_dirs[0]
        required.extend([pair / f"{source}_to_{target}.txt", pair / f"{target}_to_{source}.txt"])
    return required


def main() -> int:
    config = load_config(CONFIG)
    missing = [path.relative_to(ROOT) for path in required_inputs(config) if not path.is_file()]
    if missing:
        print("SAMap has not started because these configured inputs are missing:")
        print("\n".join(f"  - {path}" for path in missing))
        return 2

    command = [sys.executable, "-m", "workflow.scripts.run_original_samap", "--config", str(CONFIG)]
    print("Starting official SAMap for the configured species...")
    return subprocess.call(command, cwd=ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
