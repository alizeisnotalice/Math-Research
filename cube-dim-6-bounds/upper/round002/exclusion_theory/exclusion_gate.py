#!/usr/bin/env python3
"""Exact finite geometry for the cheapest body-frame cost-exclusion gate."""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from hashlib import sha256
import json
import sympy as sp


def main() -> None:
    vectors = [(1,) + s for s in product((-1, 1), repeat=5)]
    x = sp.Matrix([1, sp.Rational(3, 4), sp.Rational(3, 4), sp.Rational(3, 4), sp.Rational(3, 4), 0])
    best = None
    best_pair = None
    for i, j in combinations(range(32), 2):
        A = sp.Matrix([vectors[i], vectors[j]]).T
        residue = x - A * (A.T * A).inv() * A.T * x
        square = sp.factor((residue.T * residue)[0])
        if best is None or square < best:
            best, best_pair = square, (i, j)
    delta = sp.sqrt(best) / (60 * sp.sqrt(6))
    local_constant = sp.Rational(4960, 32) / (sp.pi**6 * delta**6)
    result = {
        "anchor": ["1", "3/4", "3/4", "3/4", "3/4", "0"],
        "distance_squared_to_S2": str(best),
        "nearest_plane_direction_indices": list(best_pair),
        "safe_delta_for_three_large_parameter_argument": str(delta),
        "generic_tail_constant_formula": str(local_constant),
        "tail_bound_at_epsilon_1_over_2000": str(sp.simplify(local_constant / 2000)),
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    assert best == sp.Rational(1, 20)
    Path(__file__).with_name("exclusion_gate.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
