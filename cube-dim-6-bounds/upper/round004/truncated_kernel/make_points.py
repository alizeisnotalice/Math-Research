#!/usr/bin/env python3
"""Freeze the three rational points for the round-004 H_epsilon pilot."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
P = 1_000_000_007
D = 4 * P
rows = [
    # The round-003 local-search minimum, scaled to common denominator D.
    (D, 4 * 998_836_847, 4 * 892_738_545, 4 * 807_464_540, 4 * 467_291_939, 4 * 96_968_636),
    # Exactly S2: absolute values are only 1 and 3/4.
    (D, 3 * P, 3 * P, 3 * P, 3 * P, 3 * P),
    # Round-002 exact Laurent/Green anchor, with a zero/repeated-node layer.
    (D, 3 * P, 3 * P, 3 * P, 3 * P, 0),
]
epsilons = [(1, 100), (1, 200), (1, 500)]


def main() -> None:
    points = HERE / "points.txt"
    points.write_text(f"{D} {len(rows)}\n" + "\n".join(" ".join(map(str, r)) for r in rows) + "\n")
    eps = HERE / "epsilons.txt"
    eps.write_text(f"{len(epsilons)}\n" + "\n".join(f"{a} {b}" for a, b in epsilons) + "\n")
    payload = {
        "status": "frozen_round004_H_epsilon_inputs",
        "scope": "3 rational representatives x 3 positive epsilons; direct truncated kernel diagnostic only",
        "points": {
            "minimum_near": list(rows[0]),
            "S2_two_level": list(rows[1]),
            "round002_green_anchor": list(rows[2]),
        },
        "epsilons": ["1/100", "1/200", "1/500"],
        "points_sha256": hashlib.sha256(points.read_bytes()).hexdigest(),
        "epsilons_sha256": hashlib.sha256(eps.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE / "input_metadata.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__": main()
