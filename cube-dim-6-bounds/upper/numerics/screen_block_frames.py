#!/usr/bin/env python3
"""Finite numerical screen for two cube-symmetric six-dimensional Cauchy frames.

This script deliberately computes *candidates*, not rigorous bounds.  It records
the quadrature choices so that any promising value can be replaced by interval
certification before it is ever used in a weak-(1,1) statement.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.stats import qmc

N = 6
Q = N / (N - 1)
SEED = 601251
ROOT = Path(__file__).resolve().parent


def log_time_rule(order: int = 96, cutoff: float = 24.0):
    """Gauss--Legendre rule for int exp(s)*p_exp(s) ds on [-cutoff,cutoff]."""
    z, w = roots_legendre(order)
    s = cutoff * z
    return np.exp(s), cutoff * w * np.exp(s)


T, TW = log_time_rule()


def density_coordinate(x: np.ndarray, t: float) -> np.ndarray:
    return (t / math.pi) ** 6 / np.prod(t * t + x * x, axis=1)


def density_three_pairs(x: np.ndarray, t: float) -> np.ndarray:
    result = np.ones(x.shape[0])
    for j in (0, 2, 4):
        a = (x[:, j] + x[:, j + 1]) / 2.0
        b = (x[:, j] - x[:, j + 1]) / 2.0
        result *= t * t / (2.0 * math.pi**2 * (t * t + a * a) * (t * t + b * b))
    return result


def green(x: np.ndarray, density) -> np.ndarray:
    value = np.zeros(x.shape[0])
    for t, wt in zip(T, TW, strict=True):
        value += wt * density(x, float(t))
    return value


def points_for_minimum() -> np.ndarray:
    # 1,024 scrambled interior points + all 64 cube corners + 729 structured
    # points.  Points on singular subspaces are retained; they can only increase
    # the sampled value and therefore never falsely improve an empirical minimum.
    interior = qmc.Sobol(d=6, scramble=True, seed=SEED).random_base2(10)
    corners = np.array(np.meshgrid(*([[0.0, 1.0]] * 6), indexing="ij")).reshape(6, -1).T
    grid = np.array(np.meshgrid(*([[0.0, 0.5, 1.0]] * 6), indexing="ij")).reshape(6, -1).T
    return np.vstack((interior, corners, grid))


def points_for_face() -> np.ndarray:
    y = qmc.Sobol(d=5, scramble=True, seed=SEED + 1).random_base2(11)
    return np.column_stack((np.ones(len(y)), y))


def screen(name: str, density) -> dict:
    xmin = points_for_minimum()
    gvals = green(xmin, density)
    finite = np.isfinite(gvals)
    gemp = float(np.min(gvals[finite]))
    argmin = xmin[np.where(gvals == gemp)[0][0]].tolist()

    xface = points_for_face()
    face_values = green(xface, density) ** Q
    Iemp = float(np.mean(face_values))
    # Cube symmetry of both selected frames gives J=6*2^6 I.  With c=g/6,
    # the formal fixed-frame cost is 6^(6/5) I / g^(6/5).
    cost = float(6.0**Q * Iemp / gemp**Q)
    return {
        "model": name,
        "frame": (
            "six coordinate Cauchy directions" if name == "coordinate" else
            "three independent (1,1),(1,-1) Cauchy direction pairs"
        ),
        "points_minimum": int(len(xmin)),
        "points_face": int(len(xface)),
        "g_empirical": gemp,
        "g_empirical_argmin": argmin,
        "I_face_q_empirical": Iemp,
        "formal_cost_empirical": cost,
        "under_64_empirically": bool(cost < 64),
        "classification": "numerical candidate only; neither g nor I has a rigorous enclosure",
    }


def source_hash() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main() -> None:
    output = {
        "purpose": "finite numerical screen; not a C6 certificate",
        "n": N,
        "q": Q,
        "quadrature": {"variable": "s=log(t)", "order": 96, "cutoff": 24.0},
        "samples": {"minimum": "1024 Sobol + 64 corners + 729 grid", "face": "2048 scrambled Sobol"},
        "cost_formula": "6^(6/5)*I/g^(6/5), valid only after cube-symmetric rigorous g and I certificates plus the full transfer proof",
        "results": [screen("coordinate", density_coordinate), screen("three_pairs", density_three_pairs)],
        "environment": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
        "source_sha256": source_hash(),
    }
    (ROOT / "screen_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
