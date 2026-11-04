#!/usr/bin/env python3
"""Numerical checks for the fixed-cube finite-scale theorem.

Sampling and quadrature check the formulas; the proof is in the paper.
Run with the project directory as cwd. No generated full search data are saved.
"""
from __future__ import annotations

import argparse
import json
import platform
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq


def cube_value(x, ell, beta):
    d = np.abs(np.asarray(x, dtype=float))
    overlap = np.maximum(0.0, np.minimum(np.minimum(1.0, ell), (1.0 + ell) / 2.0 - d))
    return float(ell ** (beta - len(d)) * np.prod(overlap))


def finite_value(x, beta):
    d = np.abs(np.asarray(x, dtype=float))
    if np.any(d >= 0.5):
        raise ValueError("The finite contact-scale theorem requires x in int(Q0).")
    scales = np.unique(np.concatenate([1.0 - 2.0 * d, 1.0 + 2.0 * d]))
    return max(cube_value(d, ell, beta) for ell in scales)


def value2(r, s, beta):
    return max((1 - 2 * r) ** beta,
               (1 - 2 * s) ** (beta - 1) * (1 - s - r),
               (1 + 2 * s) ** (beta - 2) * (1 + s - r),
               (1 + 2 * r) ** (beta - 2))


def formula2(beta, alpha):
    """Fixed-beta eventual formula, evaluated at 70 decimal digits."""
    with mp.workdps(70):
        b, a = mp.mpf(str(beta)), mp.mpf(str(alpha))
        c = 2 - b
        if b < mp.mpf("0.5"):
            v = (-mp.expm1(mp.log(a) / b)) ** 2
        elif b < 1:
            v = 2 - 4 * a / (2 - b) + 2 * b * a ** (2 / b) / (2 - b)
        elif b < mp.mpf("1.5"):
            v = 4 * (c * a ** (-1 / c) / (c + 1) - 1 + a / (c + 1))
        else:
            v = mp.expm1(-mp.log(a) / c) ** 2
        return float(v)


def area2(beta, alpha):
    """Independent boundary root integration using the four scale candidates."""
    diagonal = max(-np.expm1(np.log(alpha) / beta),
                   np.expm1(-np.log(alpha) / (2 - beta))) / 2
    def width(s):
        if s >= diagonal:
            return 0.0
        r = brentq(lambda r: value2(r, s, beta) - alpha,
                   s, 0.499999999, xtol=5e-15)
        return max(0.0, r - s)
    integral, error = quad(width, 0, diagonal, epsabs=2e-14, epsrel=2e-10, limit=150)
    return 8 * integral, 8 * error


def area_job(args):
    beta, alpha = args
    observed, error = area2(beta, alpha)
    expected = formula2(beta, alpha)
    relative = abs(observed - expected) / expected
    return {"beta": beta, "alpha": alpha, "quadrature": observed,
            "formula": expected, "quad_error_estimate": error,
            "relative_difference": relative}


def run(workers=12):
    seed = 20261104
    rng = np.random.default_rng(seed)
    max_excess = 0.0
    samples = 0
    # Interior centers, each fixed-scale value is compared with contact scales.
    for n in (1, 2, 3, 5, 8):
        for fraction in (0.05, 0.25, 0.5, 0.75, 0.95):
            beta = n * fraction
            for x in rng.uniform(-0.49, 0.49, size=(80, n)):
                maximum = finite_value(x, beta)
                for ell in np.exp(rng.uniform(np.log(0.005), np.log(5), 16)):
                    max_excess = max(max_excess, cube_value(x, ell, beta) - maximum)
                    samples += 1
    assert max_excess < 2e-12, max_excess
    betas = (0.2, 0.49, 0.5, 0.51, 0.8, 0.99, 1.0, 1.01, 1.49, 1.5, 1.51, 1.8)
    jobs = [(b, 1.0 - 2.0 ** (-j)) for b in betas for j in (8, 10, 12, 14)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(area_job, jobs))
    maximum_difference = max(row["relative_difference"] for row in rows)
    assert maximum_difference < 2e-7, maximum_difference
    return {"kind": "numerical-cross-check", "seed": seed,
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "mpmath": mp.__version__,
            "workers": workers, "contact_scale_samples": samples,
            "maximum_sampled_excess": max_excess,
            "area_cases": len(rows), "maximum_area_relative_difference": maximum_difference,
            "coverage": "Interior centers only; finite tests do not certify all inputs or thresholds.",
            "area_rows": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--output", default="results/fixed_cube_python.json")
    args = parser.parse_args()
    result = run(args.workers)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "area_rows"}, ensure_ascii=False))
