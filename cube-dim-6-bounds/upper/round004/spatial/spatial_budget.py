#!/usr/bin/env python3
"""Exact sizing table for the H_epsilon spatial-comparison certificate."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import sympy as sp


def main() -> None:
    K = sp.Rational(6, 7) * sp.Rational(32, 243) ** sp.Rational(1, 6) * sp.pi ** sp.Rational(-11, 6)
    rows = []
    for epsilon in [sp.Rational(1, 20), sp.Rational(1, 50), sp.Rational(1, 100), sp.Rational(1, 200), sp.Rational(1, 1000)]:
        cube_cells_per_axis = sp.ceiling(4 / epsilon)  # [-1,1], side 2r=epsilon/2
        face_cells_per_axis = sp.ceiling(2 / epsilon)  # [0,1], same side 2r
        rows.append({
            "epsilon": str(epsilon),
            "half_width": str(epsilon / 4),
            "tail_norm_upper": str(K * epsilon ** sp.Rational(7, 6)),
            "cube_cells": str(cube_cells_per_axis ** 6),
            "positive_face_cells": str(face_cells_per_axis ** 5),
        })
    payload = {
        "tail_coefficient": str(K),
        "comparison_radius_rule": "r=epsilon/4 gives exp(+-1/4) multiplicative H_epsilon box factors",
        "rows": rows,
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("spatial_budget.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
