#!/usr/bin/env python3
"""Reproducible small LP comparison for product rectangle catalogues.

For each base catalogue and star-product parameter set this solves (a) the
full atom LP for A x B, and (b) its symmetric layer reduction (2).  It then
converts each floating optimizer output into exact rational primal/dual
feasibility certificates and checks that the resulting objective intervals
overlap.  No large catalogue is enumerated.
"""
from __future__ import annotations

import itertools
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
TOL = 1e-9
DENOM = 10**9
MIX = F(1, 10**12)


def enc(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def rect(lo, hi):
    return [[F(a), F(b)] for a, b in zip(lo, hi)]


def base_catalogues():
    return [
        {"name": "1d_disjoint_noncommon", "dimension": 1,
         "boxes": [rect([0], [1]), rect([2], [3])],
         "class": "noncommon_disjoint"},
        {"name": "1d_overlapping_noncommon", "dimension": 1,
         "boxes": [rect([0], [2]), rect([1], [3]), rect([F(5, 2)], [4])],
         "class": "overlapping_without_common_core"},
        {"name": "1d_common_core", "dimension": 1,
         "boxes": [rect([-1], [1]), rect([-2], [1]), rect([-1], [2])],
         "class": "common_core"},
        {"name": "2d_disjoint_noncommon", "dimension": 2,
         "boxes": [rect([0, 0], [1, 1]), rect([2, 0], [3, 1])],
         "class": "noncommon_disjoint"},
        {"name": "2d_overlapping_noncommon", "dimension": 2,
         "boxes": [rect([0, 0], [2, 2]), rect([1, 0], [3, 2]),
                   rect([2, 1], [4, 3])],
         "class": "overlapping_without_common_core"},
        {"name": "2d_common_core", "dimension": 2,
         "boxes": [rect([-1, -1], [1, 1]), rect([-2, -1], [1, 1]),
                   rect([-1, -1], [2, 2])],
         "class": "common_core"},
    ]


def geometry(boxes):
    """Exact coordinate-grid atoms covering the union of a finite box family."""
    d = len(boxes[0])
    axes = [sorted({z for box in boxes for z in box[q]}) for q in range(d)]
    shape = [len(axis) - 1 for axis in axes]
    cells = []
    for idx in itertools.product(*(range(v) for v in shape)):
        lo = [axes[q][idx[q]] for q in range(d)]
        hi = [axes[q][idx[q] + 1] for q in range(d)]
        mask = [all(lo[q] >= box[q][0] and hi[q] <= box[q][1]
                    for q in range(d)) for box in boxes]
        if any(mask):
            mass = math.prod((hi[q] - lo[q] for q in range(d)), start=F(1))
            cells.append({"lo": lo, "hi": hi, "mass": mass, "mask": mask})
    volumes = [math.prod((axis[1] - axis[0] for axis in box), start=F(1))
               for box in boxes]
    for i, vol in enumerate(volumes):
        assert sum((cell["mass"] for cell in cells if cell["mask"][i]), F(0)) == vol
    return cells, volumes


def star_catalogue(k, h, a):
    """Complete h-fold product of k rows, one extended axis per row per block."""
    rows = []
    for choices in itertools.product(range(k), repeat=h):
        box = []
        for block, extended_axis in enumerate(choices):
            for axis in range(k):
                box.append([F(0), F(1) + a if axis == extended_axis else F(1)])
        rows.append(box)
    return rows


def product_boxes(base, star):
    return [[*base_box, *star_box] for base_box in base for star_box in star]


def build_full_lp(base, star, p):
    boxes = product_boxes(base, star)
    cells, volumes = geometry(boxes)
    nvar, nrow = len(cells), len(boxes)
    A = [[(cell["mass"] if cell["mask"][i] else F(0)) for cell in cells]
         for i in range(nrow)]
    b = [p * vol for vol in volumes]
    c = [cell["mass"] for cell in cells]
    return boxes, cells, A, b, c


def build_layer_lp(base_boxes, k, h, a, p):
    base_cells, base_volumes = geometry(base_boxes)
    hcoef = [F(math.comb(h, j)) * a**j for j in range(h + 1)]
    layer_cost = [hcoef[j] * k**j for j in range(h + 1)]
    v = (1 + a) ** h
    nbase, nlayers = len(base_cells), h + 1
    nvar, nrow = nbase * nlayers, len(base_boxes)
    c = [layer_cost[j] * base_cells[i]["mass"]
         for i in range(nbase) for j in range(nlayers)]
    A = []
    b = []
    for r, vol in enumerate(base_volumes):
        row = []
        for i, cell in enumerate(base_cells):
            for j in range(nlayers):
                row.append(hcoef[j] * cell["mass"] if cell["mask"][r] else F(0))
        A.append(row)
        b.append(p * v * vol)
    return base_cells, hcoef, layer_cost, v, A, b, c


def solve_lp(Aq, bq, cq, label):
    """Numerically solve, then certify with exact rational primal and dual."""
    A = np.asarray([[float(v) for v in row] for row in Aq], dtype=float)
    b = np.asarray([float(v) for v in bq], dtype=float)
    c = np.asarray([float(v) for v in cq], dtype=float)
    started = time.time()
    result = linprog(c, A_ub=-A, b_ub=-b, bounds=(0.0, 1.0), method="highs",
                     options={"primal_feasibility_tolerance": TOL,
                              "dual_feasibility_tolerance": TOL})
    elapsed = time.time() - started
    if not result.success:
        return {"label": label, "solver_status": int(result.status),
                "solver_message": str(result.message), "solver_success": False,
                "tolerance": TOL, "seconds": elapsed}

    # Exact primal: bounded-denominator rationalization, followed by a uniform
    # mixture with 1 just large enough to repair any rationalization shortfall.
    x = [min(F(1), max(F(0), F(float(v)).limit_denominator(DENOM)))
         for v in result.x]
    eps = F(0)
    for row, rhs in zip(Aq, bq):
        current = sum((coef * value for coef, value in zip(row, x)), F(0))
        capacity = sum(row, F(0))
        if current < rhs:
            assert capacity > rhs
            eps = max(eps, (rhs - current) / (capacity - current))
    if eps:
        eps += MIX
    else:
        eps = MIX
    assert eps < 1
    x = [(1 - eps) * value + eps for value in x]
    assert all(F(0) <= value <= F(1) for value in x)
    row_margins = [sum((coef * value for coef, value in zip(row, x)), F(0)) - rhs
                   for row, rhs in zip(Aq, bq)]
    assert all(margin >= 0 for margin in row_margins), f"{label}: rational primal infeasible"
    primal = sum((coef * value for coef, value in zip(cq, x)), F(0))

    # For min c.x subject to A x >= b, 0<=x<=1, any y>=0 gives lower bound
    # b.y + sum_j min(0, c_j-(A^T y)_j). Recover and rationalize HiGHS y.
    yf = [-float(value) for value in result.ineqlin.marginals]
    y = [max(F(0), F(v).limit_denominator(DENOM)) for v in yf]
    assert all(v >= 0 for v in y), f"{label}: negative dual multiplier"
    coeff = [F(v) for v in cq]
    for i, multiplier in enumerate(y):
        if multiplier:
            for j, coef in enumerate(Aq[i]):
                coeff[j] -= coef * multiplier
    dual = sum((rhs * multiplier for rhs, multiplier in zip(bq, y)), F(0))
    dual += sum((min(F(0), value) for value in coeff), F(0))
    assert dual <= primal, f"{label}: exact weak duality failed"
    return {
        "label": label,
        "solver_success": True,
        "solver_status": int(result.status),
        "solver_message": str(result.message),
        "solver_objective_float": float(result.fun),
        "tolerance": TOL,
        "seconds": elapsed,
        "variables": len(cq),
        "constraints": len(bq),
        "primal_upper": enc(primal),
        "dual_lower": enc(dual),
        "objective_interval": [enc(dual), enc(primal)],
        "rational_mixture_epsilon": enc(eps),
        "min_exact_primal_slack": enc(min(row_margins)),
        "exact_relative_gap": enc((primal - dual) / dual) if dual > 0 else None,
        "primal_rational": [enc(v) for v in x],
        "dual_rational": [enc(v) for v in y],
    }


def job(task):
    base, k, h, a, p, seed = task
    np.random.seed(seed)  # deterministic task seed; LP itself is deterministic.
    star = star_catalogue(k, h, a)
    full_boxes, full_cells, full_A, full_b, full_c = build_full_lp(base["boxes"], star, p)
    layer_cells, hcoef, layer_cost, v, layer_A, layer_b, layer_c = build_layer_lp(
        base["boxes"], k, h, a, p)
    full_result = solve_lp(full_A, full_b, full_c, "full_product")
    layer_result = solve_lp(layer_A, layer_b, layer_c, "symmetric_layers")
    ok = False
    intersect = None
    if full_result["solver_success"] and layer_result["solver_success"]:
        fi = list(map(F, full_result["objective_interval"]))
        li = list(map(F, layer_result["objective_interval"]))
        lower, upper = max(fi[0], li[0]), min(fi[1], li[1])
        ok = lower <= upper
        if ok:
            intersect = [enc(lower), enc(upper)]
    return {
        "base_name": base["name"],
        "base_class": base["class"],
        "base_dimension": base["dimension"],
        "base_boxes": [[[enc(z) for z in axis] for axis in box]
                       for box in base["boxes"]],
        "k": k, "h": h, "a": enc(a), "p": enc(p), "seed": seed,
        "star_rows": len(star),
        "full_dimension": base["dimension"] + h * k,
        "full_atoms": len(full_cells),
        "layer_atoms": len(layer_cells),
        "layer_variables": (h + 1) * len(layer_cells),
        "layer_capacities": [enc(c) for c in hcoef],
        "layer_total_volume": enc(v),
        "solver_method": "scipy.optimize.linprog(method='highs')",
        "seed_note": "fixed NumPy seed per task; LP solver is deterministic",
        "tolerance": TOL,
        "full": full_result,
        "layers": layer_result,
        "exact_intervals_overlap": ok,
        "overlap_interval": intersect,
    }


def main():
    bases = base_catalogues()
    tasks = []
    seed0 = 711202611
    for bi, base in enumerate(bases):
        for k in (2, 3):
            for h in (1, 2):
                for a in (F(1, 2), F(1), F(2)):
                    for p in (F(1, 4), F(1, 2), F(3, 4)):
                        seed = seed0 + bi * 10000 + k * 1000 + h * 100 + int(2 * a) * 10 + int(4 * p)
                        tasks.append((base, k, h, a, p, seed))
    started = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(job, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda x: (x["base_name"], x["k"], x["h"], F(x["a"]), F(x["p"])))
    successful = [r for r in results if r["full"]["solver_success"] and r["layers"]["solver_success"]]
    overlaps = [r for r in successful if r["exact_intervals_overlap"]]
    all_gaps = []
    for r in successful:
        for key in ("full", "layers"):
            val = r[key]["exact_relative_gap"]
            if val is not None:
                all_gaps.append(F(val))
    output = {
        "status": "complete" if len(successful) == len(tasks) and len(overlaps) == len(tasks) else "partial_or_failed",
        "description": "Small exact-grid base catalogues and complete k-row, h-block star products; LP optima independently enclosed by rational primal/dual certificates.",
        "case_count": len(results),
        "expected_case_count": len(tasks),
        "successful_full_and_layer_solve_count": len(successful),
        "exact_interval_overlap_count": len(overlaps),
        "each_case_has_two_exact_rational_objective_intervals": all(
            "objective_interval" in r["full"] and "objective_interval" in r["layers"] for r in successful),
        "maximum_exact_relative_gap_across_all_lp_certificates": enc(max(all_gaps, default=F(0))),
        "wall_seconds": time.time() - started,
        "workers": 12,
        "blas_threads_per_worker": 1,
        "solver_tolerances": {"primal_feasibility": TOL, "dual_feasibility": TOL},
        "numeric_solver": "SciPy HiGHS; all returned objective bounds rechecked with exact Fraction feasibility",
        "fraction_certification": "rationalized primal repaired by exact uniform mixture with 1; nonnegative rationalized dual yields exact lower bound; reported objective intervals are rigorous for each finite LP",
        "cases": results,
    }
    (HERE / "product_lp.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: output[k] for k in (
        "status", "case_count", "successful_full_and_layer_solve_count",
        "exact_interval_overlap_count", "maximum_exact_relative_gap_across_all_lp_certificates",
        "wall_seconds", "workers", "blas_threads_per_worker")}, indent=2))


if __name__ == "__main__":
    main()
