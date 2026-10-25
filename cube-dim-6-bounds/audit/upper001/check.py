#!/usr/bin/env python3
"""Independent finite algebra checks for the body-frame upper001 audit.

This checker deliberately certifies existence of the advertised determinant-160
bases, never a global maximum determinant.
"""

from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import sympy as sp


def main() -> None:
    vectors = [(1,) + signs for signs in product((-1, 1), repeat=5)]
    certified = (0, 3, 5, 14, 22, 25)
    b = sp.Matrix([vectors[i] for i in certified])
    assert abs(b.det()) == 160
    assert b.T * b != 32 * sp.eye(6)  # avoid confusing Gram of frame with basis
    frame = sp.Matrix(vectors)
    assert frame.T * frame == 32 * sp.eye(6)

    ranks = {}
    for k in (1, 2, 3):
        values = [sp.Matrix([vectors[i] for i in ids]).rank()
                  for ids in combinations(range(32), k)]
        assert min(values) == max(values) == k
        ranks[k] = len(values)

    # The diagonal group flips coordinates 2,...,6.  It is transitive on V,
    # preserves |det|, and maps the certified basis (which contains vector 0)
    # to a determinant-160 basis containing every desired active direction.
    completions = []
    for target in vectors:
        diagonal = sp.diag(1, *[target[j] * vectors[0][j] for j in range(1, 6)])
        transformed = [tuple(int(x) for x in (diagonal * sp.Matrix(row))) for row in (vectors[i] for i in certified)]
        assert target in transformed
        det = abs(sp.Matrix(transformed).det())
        assert det == 160
        completions.append(det)

    pi = sp.pi
    cf = sp.Rational(32)**3 / (160 * pi**6)
    assert sp.simplify(cf - sp.Rational(1024, 5) / pi**6) == 0
    # Eq. (14): 6/7 * 32^(-5/6) * Cf^(1/6) * (32/pi)^(5/6).
    lhs = sp.Rational(6, 7) * 32**sp.Rational(-5, 6) * cf**sp.Rational(1, 6) * (32/pi)**sp.Rational(5, 6)
    rhs = 6 * 2**sp.Rational(5, 3) / (7 * 5**sp.Rational(1, 6) * pi**sp.Rational(11, 6))
    assert sp.simplify(lhs / rhs - 1) == 0
    # Exact cancellation in the conditional cost calculation.
    g, I = sp.symbols('g I', positive=True)
    c = g / 6
    J = 6 * 64 * I
    mass = sp.Rational(5, 6) * c**sp.Rational(-1, 5) * J / (64 * (g - c))
    assert sp.simplify(mass - 6**sp.Rational(6, 5) * I / g**sp.Rational(6, 5)) == 0

    payload = {
        "status": "PASS",
        "exact_det160_basis": list(certified),
        "subset_rank_counts": ranks,
        "all_active_directions_have_det160_completion": len(completions),
        "Cf": "1024/(5*pi^6)",
        "whole_face_epsilon_coefficient": "6*2^(5/3)/(7*5^(1/6)*pi^(11/6))",
        "conditional_cost": "6^(6/5)*I/g0^(6/5)",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
