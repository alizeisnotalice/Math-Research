"""Round 8 finite numerical study of psi_2 and psi_3 in the Cauchy kernel.

Closed slice formulas come from 2 PV integral f'(y)/(y-u) for
f(y)=(a-|y|)_+^p.  The LP enforces the true two-variable generator formulas
on a finite r-by-angle grid and radial kernel domination on the same r grid.
It is discovery only, never a general weak-type certificate.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parent
R = 4746319 / 2500000
INCUMBENT = 27226955153507 / 7019746875000
N_INNER, N_OUTER, NA = 100, 50, 30  # 4,500 joint points, plus 100 dominance points
B = 2.0 - 2.0 * math.log(2.0)


def h2(x: float) -> float:
    """jumpGen1 of (1-|x|)_+^2 for x>0 (including the exterior x>1)."""
    if x == 1.0:
        return 4 * (2 - 2 * math.log(2))
    return 4 * (2 - (1 - x) * math.log(abs(1 - x))
                - (1 + x) * math.log(1 + x) + 2 * math.log(x))


def h3(x: float) -> float:
    """jumpGen1 of (1-|x|)_+^3 for x>0, from polynomial division plus PV logs."""
    if x == 1.0:
        return 6 * (3 - 4 * math.log(2))
    return 6 * (3 + (1 - x) ** 2 * math.log(x / abs(1 - x))
              + (1 + x) ** 2 * math.log(x / (1 + x)))


def slice_generator_numeric(p: int, x: float) -> float:
    """Direct finite-interval integral of the original second-difference formula."""
    def f(y: float) -> float:
        return max(1 - abs(y), 0.0) ** p
    def integrand(t: float) -> float:
        if abs(t) < 1e-7:
            return p * (p - 1) * (1.0 - x) ** (p - 2) if x < 1.0 else 0.0
        return (f(x + t) + f(x - t) - 2.0 * f(x)) / (t * t)
    # The integrand vanishes outside [-(1+x), 1+x]; pass all kinks to quad.
    end = 1 + x
    kinks = sorted({-x - 1, -x, -x + 1, x - 1, x, x + 1, 0.0})
    core = quad(integrand, -end, end, points=[t for t in kinks if -end < t < end],
                epsabs=2e-10, epsrel=2e-10, limit=300)[0]
    # Beyond end both translated terms vanish, but the -2 f(x) term remains.
    return core - 4.0 * f(x) / end


def tail_f(c: float) -> float:
    if c == 0.0:
        return 1.0 / (2.0 * R)
    return -(R / c**2) * math.log1p(-c / R) - 1.0 / c


def tail_generator(r: float, d: float) -> float:
    return 2.0 * (2.0 * tail_f(r) + tail_f(d) + tail_f(-d)) / (R - 1.0)


def bump_generator(p: int, u: float, v: float) -> float:
    """Full 2D jump generator of psi_p at u,v>0; includes r>=1 exterior slices."""
    h = h2 if p == 2 else h3
    def term(a: float, x: float) -> float:
        return 0.0 if a <= 0.0 else a ** (p - 1) * h(x / a)
    return term(1.0 - v, u) + term(1.0 - u, v)


def cusp(r: float) -> float:
    return (R / r - 1.0) / (R - 1.0)


def run_lp() -> dict:
    inner_radii = np.linspace(0.005, 0.995, N_INNER)
    outer_radii = np.linspace(1.002, R - 0.002, N_OUTER)
    radii = np.r_[inner_radii, outer_radii]
    angles = np.arange(NA, dtype=float) / NA  # d/r in [0,1), includes the diagonal
    rows, rhs, meta = [], [], []
    for r in radii:
        for q in angles:
            u, v = r * (1 + q) / 2, r * (1 - q) / 2
            rows.append([bump_generator(2, u, v), bump_generator(3, u, v)])
            rhs.append(tail_generator(r, r * q))
            meta.append({"kind": "joint_generator", "r": float(r), "d_over_r": float(q)})
    # Exact radial kernel comparison, sampled only for this numerical screen.
    for r in inner_radii:
        rows.append([(1 - r) ** 2, (1 - r) ** 3])
        rhs.append(cusp(r) - 1.0)
        meta.append({"kind": "dominance", "r": float(r), "d_over_r": None})
    # This is a rigorous all-angle exterior sufficient condition, not a sample:
    # for r >= 1, psi_p(z)=0 and 0 <= psi_3 <= psi_2, hence
    # A psi_p <= 4 b, while A H_R >= 4/[R(R-1)].
    rows.append([1.0, 1.0])
    rhs.append(1.0 / (R * (R - 1.0) * B))
    meta.append({"kind": "exterior_global", "r": None, "d_over_r": None})
    A, b = np.asarray(rows), np.asarray(rhs)
    result = linprog(c=[-1 / 6, -1 / 10], A_ub=A, b_ub=b,
                     bounds=[(0, None), (0, None)], method="highs")
    if not result.success:
        raise RuntimeError(result.message)
    e2, e3 = map(float, result.x)
    joint_count = len(radii) * NA
    generator_slack = np.asarray(rhs[:joint_count]) - np.asarray(rows[:joint_count]) @ result.x
    dominance_slack = np.asarray(rhs[joint_count:]) - np.asarray(rows[joint_count:]) @ result.x
    dual = -np.asarray(result.ineqlin.marginals)  # convert minimization dual to max-benefit sign
    active = sorted(({
        **meta[i], "slack": float(np.asarray(rhs)[i] - np.asarray(rows)[i] @ result.x),
        "maximization_dual": float(dual[i]),
    } for i in range(len(rhs))), key=lambda row: (row["slack"], -row["maximization_dual"]))[:12]
    cost = R * R / (R - 1) - e2 / 6 - e3 / 10
    return {
        "R": R, "joint_generator_grid": {"inner": [N_INNER, NA], "exterior": [N_OUTER, NA]},
        "joint_generator_points": joint_count, "dominance_points": N_INNER,
        "total_linear_constraints": len(rhs),
        "e2": e2, "e3": e3, "cost_candidate": cost,
        "incumbent": INCUMBENT, "candidate_minus_incumbent": cost - INCUMBENT,
        "min_joint_generator_slack": float(generator_slack.min()),
        "min_dominance_slack": float(dominance_slack.min()),
        "most_active_constraints": active,
        "linprog_status": result.message,
        "status": "finite_grid_numerical_candidate_only",
    }


def main() -> None:
    cross = []
    for p, formula in ((2, h2), (3, h3)):
        for x in (0.17, 0.43, 0.79, 1.23):
            numerical = slice_generator_numeric(p, x)
            cross.append({"p": p, "x": x, "closed": formula(x), "direct_integral": numerical,
                          "absolute_error": abs(formula(x) - numerical)})
    worst_cross_error = max(row["absolute_error"] for row in cross)
    if worst_cross_error > 5e-7:
        raise AssertionError(f"slice formula/direct-integral discrepancy: {worst_cross_error}")
    out = {"slice_cross_checks": cross, "worst_slice_cross_error": worst_cross_error, "lp": run_lp()}
    (ROOT / "search_result.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
