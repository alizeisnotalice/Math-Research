#!/usr/bin/env python3
"""Round 14: fixed 6000-point interior LPs with the proved exterior row.

Numerical discovery only. At most 18 LPs; duplicate rational R values reuse
their existing solution. Convex bump values are evaluated once and reused.
No continuous certificate is attempted by this program.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.optimize import linprog

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / "round012" / "search_quartic.py"
spec = importlib.util.spec_from_file_location("round014_round012_numeric_core", SOURCE)
core = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = core
spec.loader.exec_module(core)
INCUMBENT = Q(173331, 44750)
LP_BUDGET = 18
POINT_BUDGET = 6000
COUNTS = {"distinct_interior_points": 0, "Apsi2": 0, "Apsi3": 0,
          "Apsi4": 0, "AH_R": 0, "LPs": 0}
H = np.array([8-8*math.log(2), 18-24*math.log(2), 136/3-64*math.log(2)])
C = np.array([-1/6, -1/10, -1/15])


def fixed_points():
    regular = np.linspace(.0025, .9975, 180)
    extra = np.array(list(map(float, [
        ".000001", ".00001", ".0001", ".0005", ".001", ".80", ".82",
        ".84", ".85", ".86", ".88", ".9", ".998", ".999", ".9995",
        ".9999", ".99999", ".999999", ".99999999", ".9999999999"])))
    radii = np.sort(np.r_[regular, extra])
    angles = np.r_[np.arange(25)/25., .98, .99, .999, .99999, 1-1e-8]
    points = [(float(r), float(q)) for r in radii for q in angles]
    assert len(points) == len(set(points)) == POINT_BUDGET
    assert all(0 < r < 1 and 0 <= q < 1 for r, q in points)
    return points, radii, angles


def tail_f(R, c):
    # Round012 formula and stability branch, with R made an explicit input.
    y = c/R
    if abs(y) < .25:
        total = 0.
        power = 1.
        for n in range(40):
            total += power/((n+2)*R)
            power *= y
        return total
    return -(R/c**2)*math.log1p(-c/R)-1/c


def tail_generator(R, r, d):
    return 2*(2*tail_f(R, r)+tail_f(R, d)+tail_f(R, -d))/(R-1)


def description(i, points):
    if i == len(points):
        return {"kind": "proved_global_exterior_sufficient"}
    if i == len(points)+1:
        return {"kind": "global_kernel_dominance_sufficient"}
    r, q = points[i]
    return {"kind": "fixed_interior_generator", "r": r, "d_over_r": q,
            "u": r*(1+q)/2, "v": r*(1-q)/2}


def solve(Rq, rows, points, stage):
    R = float(Rq)
    rhs = np.array([tail_generator(R, r, r*q) for r, q in points])
    COUNTS["AH_R"] += len(points)
    D1 = 2*(2*tail_f(R, 1)+1/R)/(R-1)
    extra = np.vstack([H, np.full(3, (R-1)/4)])
    extra_rhs = np.array([D1, R])
    A = np.vstack([rows, extra])
    b = np.r_[rhs, extra_rhs]
    if COUNTS["LPs"] >= LP_BUDGET:
        raise RuntimeError("LP budget exceeded")
    COUNTS["LPs"] += 1
    result = linprog(C, A_ub=A, b_ub=b, bounds=[(0, None)]*3, method="highs",
                     options={"primal_feasibility_tolerance": 1e-9,
                              "dual_feasibility_tolerance": 1e-9})
    out = {"R_exact": str(Rq), "R": R, "stage_first_solved": stage,
           "success": bool(result.success), "solver_status": int(result.status),
           "solver_message": result.message, "iterations": int(result.nit),
           "fixed_points": len(points), "constraints": len(b),
           "D_R1": D1}
    if not result.success:
        return out
    x = result.x
    slack = b-A@x
    dual = np.asarray(result.ineqlin.marginals)
    selected = set(map(int, np.argsort(slack)[:12]))
    selected.update(map(int, np.flatnonzero(np.abs(dual) > 1e-10)))
    active = [{**description(i, points), "slack": float(slack[i]),
               "benefit_dual": float(-dual[i]), "coefficients": A[i].tolist(),
               "rhs": float(b[i])} for i in sorted(selected, key=lambda i: slack[i])]
    cost = R*R/(R-1)+float(C@x)
    out.update({"e2": float(x[0]), "e3": float(x[1]), "e4": float(x[2]),
                "cost_candidate": cost,
                "candidate_minus_strict_incumbent": cost-float(INCUMBENT),
                "minimum_grid_generator_slack": float(np.min(slack[:len(points)])),
                "global_exterior_slack": float(slack[-2]),
                "kernel_dominance_slack": float(slack[-1]),
                "positive_dual_constraints": [v for v in active if v["benefit_dual"] > 1e-10],
                "near_active_constraints": active,
                "benefit_primal": float(-result.fun),
                "benefit_dual": float(-b@dual),
                "numeric_dual_gap": float(result.fun-b@dual),
                "lower_bound_reduced_costs": np.asarray(result.lower.marginals).tolist()})
    return out


def save(data):
    data["actual_budget"] = dict(COUNTS)
    data["actual_budget"]["generator_evaluations_total"] = sum(COUNTS[k] for k in ("Apsi2", "Apsi3", "Apsi4", "AH_R"))
    (ROOT / "search_result.json").write_text(json.dumps(data, indent=2, ensure_ascii=False)+"\n")


def main():
    if (ROOT / "search_result.json").exists():
        raise SystemExit("Round14 search_result.json already exists; refusing duplicate budget use")
    started = time.monotonic()
    points, radii, angles = fixed_points()
    COUNTS["distinct_interior_points"] = len(points)
    rows = np.array([[core.bump_generator(p, r*(1+q)/2, r*(1-q)/2)
                      for p in (2, 3, 4)] for r, q in points])
    for key in ("Apsi2", "Apsi3", "Apsi4"):
        COUNTS[key] = len(points)
    if not np.isfinite(rows).all():
        raise ValueError("Nonfinite bump generator coefficient")
    arguments = [x/a for r, q in points
                 for a, x in ((1-r*(1-q)/2, r*(1+q)/2),
                              (1-r*(1+q)/2, r*(1-q)/2))]
    data = {"round": 14, "status": "running_finite_grid_only",
            "strict_incumbent": str(INCUMBENT), "LP_budget": LP_BUDGET,
            "distinct_interior_point_budget": POINT_BUDGET,
            "family": "R plus e2,e3,e4>=0; fixed 6000 interior constraints; proved exterior sufficient row and kernel dominance row",
            "source": {"path": str(SOURCE), "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                       "reused": "hp and bump_generator imported unchanged from round012; tail_f stability branch reproduced with explicit R"},
            "grid": {"radii": radii.tolist(), "d_over_r": angles.tolist(),
                     "point_digest_sha256": hashlib.sha256(json.dumps(points).encode()).hexdigest(),
                     "min_slice_argument": min(arguments), "max_slice_argument": max(arguments),
                     "all_slice_arguments_strictly_between_zero_and_one": all(0<t<1 for t in arguments)},
            "H2_H3_H4": H.tolist(),
            "environment": {"python": sys.version, "platform": platform.platform(),
                            "numpy": np.__version__, "scipy": scipy.__version__},
            "LPs": [], "R_requests": []}
    cache = {}

    def request(Rq, stage):
        key = str(Rq)
        reused = key in cache
        data["R_requests"].append({"R_exact": key, "stage": stage, "reused": reused})
        if not reused:
            entry = solve(Rq, rows, points, stage)
            cache[key] = entry
            data["LPs"].append(entry)
            save(data)
            print(json.dumps({k: entry.get(k) for k in
                  ("R_exact", "success", "e2", "e3", "e4", "cost_candidate")}), flush=True)
        return cache[key]

    stage1_r = [Q(9, 5)+Q(i, 40) for i in range(9)]
    stage1 = [request(Rq, 1) for Rq in stage1_r]
    feasible = [entry for entry in stage1 if entry["success"]]
    if not feasible:
        data["status"] = "all_first_stage_LPs_failed_no_expansion"
        save(data)
        return
    first_best = min(feasible, key=lambda v: v["cost_candidate"])
    center = Q(first_best["R_exact"])
    stage2_r = [center+Q(i-4, 160) for i in range(9)]
    stage2 = [request(Rq, 2) for Rq in stage2_r]
    best = min((entry for entry in data["LPs"] if entry["success"]), key=lambda v: v["cost_candidate"])
    data.update({"status": "finite_numerical_stage_complete_no_certificate",
                 "first_stage_best_R": str(center),
                 "first_stage_best_at_R_boundary": center in (stage1_r[0], stage1_r[-1]),
                 "second_stage_R_range": [str(stage2_r[0]), str(stage2_r[-1])],
                 "best_R_at_second_stage_boundary": Q(best["R_exact"]) in (stage2_r[0], stage2_r[-1]),
                 "best_R_at_overall_scanned_boundary": best["R"] in (min(v["R"] for v in data["LPs"]), max(v["R"] for v in data["LPs"])),
                 "best_candidate": best,
                 "elapsed_seconds": time.monotonic()-started,
                 "scope": "No strict certificate, no rational candidate selection, no new upper bound; stop after prescribed finite R scan"})
    save(data)
    print(json.dumps({"status": data["status"], "best_R": best["R_exact"],
                      "best_cost": best["cost_candidate"], "actual_budget": data["actual_budget"],
                      "elapsed_seconds": data["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
