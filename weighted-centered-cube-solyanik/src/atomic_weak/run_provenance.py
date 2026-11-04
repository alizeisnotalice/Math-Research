"""Versioned checkpoints for the round-25 exploratory programs."""
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy


CACHE_SCHEMA = "round25-signed-cache-v1"


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signature(script_path, task, model, imported_sources=()):
    """Return a deterministic signature for source, dependencies and parameters."""
    payload = {
        "cache_schema": CACHE_SCHEMA,
        "script_sha256": sha256_file(script_path),
        "imported_source_sha256": {str(Path(path).name): sha256_file(path)
                                   for path in imported_sources},
        "task": task,
        "model": model,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return {"sha256": hashlib.sha256(encoded).hexdigest(), "payload": payload}


def load_matching_or_quarantine(output, expected):
    """Return a matching record; preserve rather than silently reuse mismatches."""
    output = Path(output)
    if not output.exists():
        return None
    try:
        record = json.loads(output.read_text())
    except (OSError, json.JSONDecodeError):
        record = None
    if record and record.get("run_signature", {}).get("sha256") == expected["sha256"]:
        return record
    stale = output.parent / "stale"
    stale.mkdir(exist_ok=True)
    target = stale / (output.stem + ".signature-mismatch.json")
    suffix = 1
    while target.exists():
        target = stale / (output.stem + f".signature-mismatch-{suffix}.json")
        suffix += 1
    output.replace(target)
    return None


def write_record(output, record, run_signature):
    output = Path(output)
    record = dict(record)
    record["run_signature"] = run_signature
    temporary = output.with_suffix(".tmp")
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False))
    temporary.replace(output)
    return record
