#!/usr/bin/env python3
"""Variance-stop pilot for the 32-body-direction six-dimensional Cauchy frame.

The estimator is intentionally not a certificate.  It uses a determinant-160
six-direction base and conditions on the other 26 Cauchy coordinates.  Two
independent 512-point scrambled Sobol runs share the same finite spatial test
set; their disagreement is the stopping diagnostic.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from itertools import product
from pathlib import Path

import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parent
N = 6
Q = 6 / 5
SAMPLES_PER_REPLICATE = 512
TIME_ORDER = 80
TIME_CUTOFF = 18.0


def all_directions() -> np.ndarray:
    return np.asarray([(1,) + signs for signs in product((-1, 1), repeat=5)], dtype=float)


# Columns form an integer base with determinant -160.  The ordered row indices
# in all_directions() are (0,3,5,14,22,25), verified exactly below.
BASE_INDICES = (0, 3, 5, 14, 22, 25)


def exact_det(a: np.ndarray) -> int:
    # Bareiss elimination is sufficient for this six-by-six integer matrix.
    a = [[int(x) for x in row] for row in a]
    sign, prev = 1, 1
    for k in range(len(a) - 1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        p = a[k][k]
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                a[i][j] = (a[i][j] * p - a[i][k] * a[k][j]) // prev
        prev = p
        for i in range(k + 1, len(a)):
            a[i][k] = 0
    return sign * a[-1][-1]


V = all_directions()
B = V[list(BASE_INDICES)].T
DET_B = abs(exact_det(B.T))
assert DET_B == 160
BINV = np.linalg.inv(B)
REMAINDER = np.delete(V, BASE_INDICES, axis=0)


def time_rule():
    z, w = roots_legendre(TIME_ORDER)
    s = TIME_CUTOFF * z
    return np.exp(s), TIME_CUTOFF * w * np.exp(s)


T, TW = time_rule()


def cauchy_sum(seed: int) -> np.ndarray:
    u = qmc.Sobol(d=26, scramble=True, seed=seed).random_base2(9)
    # Scrambled Sobol never uses endpoints; tan is the exact inverse Cauchy CDF.
    c = np.tan(math.pi * (u - 0.5))
    return c @ REMAINDER


def spatial_points() -> tuple[np.ndarray, np.ndarray]:
    face = np.column_stack((np.ones(32), qmc.Sobol(d=5, scramble=True, seed=601253).random_base2(5)))
    orbits = np.array([[1.0] * k + [0.0] * (6 - k) for k in range(7)])
    # Generate a full base-2 net, then retain a fixed 25-point diagnostic
    # subsample; this avoids asking SciPy for a non-power-of-two Sobol draw.
    inside = qmc.Sobol(d=6, scramble=True, seed=601254).random_base2(5)[:25]
    cube = np.vstack((orbits, inside))
    assert len(face) == 32 and len(cube) == 32
    return face, cube


def p_base_1(z: np.ndarray) -> np.ndarray:
    """Density of B times six independent standard Cauchy variables."""
    coords = z @ BINV.T
    return np.prod(1.0 / (math.pi * (1.0 + coords * coords)), axis=1) / DET_B


def green_conditional(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Estimate E_Y int t^-6 p_B,1(x/t-Y) dt at all points x."""
    answer = np.zeros(len(x))
    for t, wt in zip(T, TW, strict=True):
        # Keep this finite array small: 32 spatial by 512 conditional samples.
        z = x[:, None, :] / t - y[None, :, :]
        base = p_base_1(z.reshape(-1, N)).reshape(len(x), len(y))
        answer += wt * (t ** -6) * np.mean(base, axis=1)
    return answer


def source_hash() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main() -> None:
    face, cube = spatial_points()
    replicas = []
    for seed in (601255, 601256):
        y = cauchy_sum(seed)
        gf = green_conditional(face, y)
        gc = green_conditional(cube, y)
        g = float(np.min(np.concatenate((gf, gc))))
        i = float(np.mean(gf**Q))
        replicas.append({
            "sobol_seed": seed,
            "g_sample_min": g,
            "I_face_sample": i,
            "formal_cost_sample": float(6.0**Q * i / g**Q),
        })
    costs = np.array([r["formal_cost_sample"] for r in replicas])
    out = {
        "purpose": "conditional Monte Carlo variance-stop pilot; not a C6 certificate",
        "frame": "32 body directions (1,eps2,...,eps6)",
        "base_indices": BASE_INDICES,
        "base_columns": B.astype(int).tolist(),
        "determinant_exact": DET_B,
        "conditional_identity": "G(x)=E_Y int_0^infty t^-6 p_B,1(x/t-Y) dt; Y=sum of other 26 directions times standard Cauchy variables",
        "analytic_checks": {
            "base_density_normalization": "exactly one by invertible linear change of variables",
            "base_density_at_zero": f"1/(160*pi^6) = {1/(160*math.pi**6):.17g}",
            "scaling": "p_B,t(x)=t^-6 p_B,1(x/t), used in the estimator",
        },
        "finite_batch": {"spatial": "32 face + 32 cube points", "conditional_samples": "two independent 512-point scrambled Sobol replicates", "time": {"log_rule_order": TIME_ORDER, "cutoff": TIME_CUTOFF}},
        "replicates": replicas,
        "cost_range": [float(costs.min()), float(costs.max())],
        "relative_spread": float((costs.max() - costs.min()) / costs.mean()),
        "stop_rule": "if replica spread is material, do not expand this estimator; derive an all-cone or variance-controlled representation first",
        "source_sha256": source_hash(),
        "environment": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__},
    }
    (ROOT / "body32_conditional_pilot_results.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
