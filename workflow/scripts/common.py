"""Shared configuration and provenance helpers."""
import hashlib, json
from pathlib import Path
import yaml

def load_config(path):
    with Path(path).open(encoding="utf-8") as handle: return yaml.safe_load(handle)

def sha256(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""): digest.update(chunk)
    return digest.hexdigest()

def write_json(path, payload):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

def fasta_ids(path):
    with Path(path).open(encoding="utf-8") as handle:
        return {line[1:].strip().split()[0] for line in handle if line.startswith(">")}

