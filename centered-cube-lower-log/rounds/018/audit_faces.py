#!/usr/bin/env python3
"""Enumerate all 2^3 supports for the specified nonnegative quadratic program.

NumPy float64 computations; this audits the supplied finite-grid covariance,
not an exact continuum covariance or an MGF remainder estimate.
"""
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "017" / "results.json"
TOL = 1e-12


def audit(record):
    C = np.asarray(record["covariance"], dtype=float)
    b = np.diag(C).copy()
    assert C.shape == (3, 3)
    assert np.max(np.abs(C - C.T)) < TOL
    eigenvalues = np.linalg.eigvalsh(C)
    assert eigenvalues[0] > 0
    faces = []
    for bits in itertools.product((0, 1), repeat=3):
        support = np.flatnonzero(bits)
        inactive = np.flatnonzero(1 - np.asarray(bits))
        theta = np.zeros(3)
        if len(support):
            theta[support] = np.linalg.solve(C[np.ix_(support, support)], b[support])
        gradient = C @ theta - b
        feasible = bool(np.all(theta >= -TOL))
        free_residual = float(np.max(np.abs(gradient[support]), initial=0))
        inactive_violation = float(np.max(-gradient[inactive], initial=0))
        kkt = feasible and free_residual < TOL and inactive_violation < TOL
        faces.append(dict(
            support_1based=(support + 1).tolist(), theta=theta.tolist(),
            objective=float(theta @ C @ theta / 2 - b @ theta),
            gradient=gradient.tolist(), primal_feasible=feasible,
            free_gradient_residual=free_residual,
            inactive_dual_violation=inactive_violation,
            complementarity_residual=float(np.max(np.abs(theta * gradient))),
            kkt=kkt))
    winners = [f for f in faces if f["kkt"]]
    assert len(winners) == 1
    best = winners[0]
    assert best["objective"] <= min(f["objective"] for f in faces if f["primal_feasible"]) + TOL
    theta = np.asarray(best["theta"])
    x = C @ theta
    primal_rate = float(x @ np.linalg.solve(C, x) / 2)
    assert np.min(x - b) >= -TOL
    assert abs(primal_rate + best["objective"]) < TOL
    unconstrained = np.linalg.solve(C, b)
    return dict(
        name=record["name"], G=record["G"], a=record["a"],
        covariance=C.tolist(), threshold_b=b.tolist(),
        eigenvalues=eigenvalues.tolist(),
        unconstrained_theta=unconstrained.tolist(),
        unconstrained_objective=float(-b @ unconstrained / 2),
        faces=faces, optimum=best,
        constrained_rate=-best["objective"],
        primal_optimum_x=x.tolist(), primal_rate=primal_rate,
        normalized_rate_to_sum_marginal_rates=float(-best["objective"] / (b.sum() / 2)),
        marginal_rate_sum=float(b.sum() / 2))


records = json.loads(SOURCE.read_text())
selected = [r for r in records if r["G"] == 131072 and r["a"] == 0.0001]
assert sorted(r["name"] for r in selected) == ["far", "near"]
result = dict(source=str(SOURCE), numpy_version=np.__version__, tolerance=TOL,
              cases=[audit(r) for r in selected])
(HERE / "audit_results.json").write_text(json.dumps(result, indent=2) + "\n")
for case in result["cases"]:
    print(case["name"])
    for key in ("unconstrained_theta", "unconstrained_objective", "optimum", "constrained_rate", "primal_optimum_x", "normalized_rate_to_sum_marginal_rates"):
        print(f"{key}: {case[key]}")
