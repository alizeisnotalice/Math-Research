#!/usr/bin/env python3
"""Freeze the five symmetry representatives for the first local-box pilot."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
D, HN = 1000, 1
rows = []
for k in range(5):
    rows.append((D, *((751,) * k), *((749,) * (4-k)), HN))


def main() -> None:
    path = HERE / "box_centers.txt"
    path.write_text(f"{D} {len(rows)}\n" + "\n".join(" ".join(map(str, row)) for row in rows) + "\n")
    payload = {
        "status": "frozen_local_box_centers",
        "center_family": "(1,(3/4+1/1000)^k,(3/4-1/1000)^(4-k),1/1000), k=0..4",
        "box": "[3/4-1/1000,3/4+1/1000]^4 x [0,1/1000] on x1=1 face",
        "epsilon": "1/100",
        "denominator": D,
        "rows": [list(row) for row in rows],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "centers_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    (HERE / "input_metadata.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__": main()
