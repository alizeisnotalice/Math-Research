#!/usr/bin/env python3
"""n=2,3 local mass checks for norm-infinity radial singularities.

For w_j(x)=c_j ||x-z_j||_infty^(tau_j-n), the mass of the centred cube of
radius r is exactly c_j*2^n*n*r^tau/tau.  Other mixture terms are evaluated
by tensor Gauss-Legendre quadrature on a fixed rational cube.  This checks the
expected O(epsilon^(n-tau_j)) local relative error when tau_j<n.
"""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parent


def exact_norminf_mass(n, c, tau, radius):
    return c * (2**n) * n * radius**tau / tau


def tensor_integral(n, component, center, radius, order=24):
    kind, c, z, tau = component
    nodes, weights = leggauss(order)
    grids = np.meshgrid(*([nodes] * n), indexing="ij")
    wg = np.meshgrid(*([weights] * n), indexing="ij")
    weight = np.prod(np.stack(wg), axis=0)
    if kind == "base":
        value = c
    else:
        coords = np.stack([center[q] + radius*grids[q] - z[q] for q in range(n)], axis=0)
        value = c * np.max(np.abs(coords), axis=0)**(tau-n)
    return float(np.sum(weight * value) * radius**n)


def run_case(n, tau):
    z0 = np.zeros(n)
    z1 = np.array([1.3] + [0.0]*(n-1))
    components = [("base", 0.35, None, None), ("singular", 0.8, z0, tau),
                  ("singular", 0.55, z1, 0.75)]
    rows = []
    # Both fixed rational scaled boxes are evaluated: E radius epsilon/2 and
    # H radius epsilon.  The local comparison concerns their individual masses.
    for k in range(3, 10):
        eps = 2.0**(-k)
        boxes = {}
        for label, radius in (("E", eps/2), ("H", eps)):
            pure = exact_norminf_mass(n, .8, tau, radius)
            approximations = [tensor_integral(n, comp, z0, radius, 12) for comp in (components[0], components[2])]
            controls = [tensor_integral(n, comp, z0, radius, 24) for comp in (components[0], components[2])]
            controls48 = [tensor_integral(n, comp, z0, radius, 48) for comp in (components[0], components[2])]
            mixed = pure + sum(controls)
            rel = abs(mixed/pure - 1)
            boxes[label] = {"radius":radius,"pure_exact":pure,"mixed_gauss24":mixed,
                "relative_error":rel,"error_over_epsilon_power":rel/eps**(n-tau),
                "gauss12_to24_difference":abs(sum(approximations)-sum(controls)),
                "gauss24_to48_difference":abs(sum(controls)-sum(controls48))}
        rows.append({"k": k, "epsilon": eps, "boxes":boxes})
    return {"n": n, "tau": tau, "condition_tau_lt_n": tau < n,
            "box_E_radius": "epsilon/2", "box_H_radius": "epsilon",
            "quadrature": "tensor Gauss-Legendre order 24 for nonsingular terms",
            "rows": rows}


def main():
    results = [run_case(n, tau) for n in (2,3) for tau in (.15,.25,.5,.75,.9)]
    output = {"scope":"local norm-infinity mass only; no nD halo optimisation",
              "results":results}
    (ROOT / "nd_blowup_summary.json").write_text(json.dumps(output, indent=2, allow_nan=False))
    for item in results:
        last = item["rows"][-1]["boxes"]["H"]
        print(f"n={item['n']} tau={item['tau']} condition={item['condition_tau_lt_n']} "
              f"last_rel={last['relative_error']:.3e} scaled={last['error_over_epsilon_power']:.3e}")


if __name__ == "__main__":
    main()
