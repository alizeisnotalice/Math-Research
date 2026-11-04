#!/usr/bin/env python3
"""Exact algebra checks for the audited original results, not proof substitutes."""
import json
import platform
from pathlib import Path

import sympy as s


def run():
    r, a, theta, t, beta, n = s.symbols("r a theta t beta n", positive=True)
    residual = s.simplify(s.diff(1 + r - 2 * a * r ** (1-theta), r, 2)
                          - 2*a*theta*(1-theta)*r**(-1-theta))
    assert residual == 0
    rect_first = s.simplify(s.diff(1 + 2*n*((1-t)**(-n/(n-beta))-1), t).subs(t, 0))
    assert s.simplify(rect_first - 2*n**2/(n-beta)) == 0
    checked_volumes = 0
    for dimension in range(1, 10):
        for sigma in [s.Rational(j, 4) for j in range(1, 4*dimension+1)]:
            k = int(s.floor(sigma))
            if k == dimension:
                continue
            # Ordered chamber integral via the Euler beta identity.
            integral = (s.factorial(dimension-k-1)*s.factorial(k)
                        / s.factorial(dimension) / sigma**(dimension-k))
            result = 2**dimension*s.factorial(dimension)*integral / (
                s.factorial(dimension-k-1)*s.factorial(k)**2)
            target = s.Rational(2**dimension, s.factorial(k))/sigma**(dimension-k)
            assert s.simplify(result-target) == 0
            checked_volumes += 1
    # Fixed-scale covariogram of the unit square; length scale is t, not sqrt(t).
    covariogram_area = 4*(t+(1-t)*s.log(1-t))
    expansion = s.series(covariogram_area, t, 0, 5)
    assert s.simplify(s.limit(covariogram_area/t**2, t, 0)-2) == 0
    return {"kind": "exact-algebra-check", "python": platform.python_version(),
            "sympy": s.__version__, "density_convexity_residual": str(residual),
            "rectangle_first_coefficient": str(rect_first),
            "tangent_volume_cases": checked_volumes,
            "square_covariogram_series": str(expansion),
            "scope": "Algebraic identities under manuscript assumptions; geometric arguments are in the paper."}


if __name__ == "__main__":
    result = run()
    target = Path("results/baseline.json")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False))
