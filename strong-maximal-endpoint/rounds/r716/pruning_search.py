#!/usr/bin/env python3
"""Finite random search: exhaustive subcatalogue pruning vs responsibility price.

This is a bounded diagnostic search, not a proof of a dimension-uniform claim.
Each catalogue uses 6 or 7 integer-endpoint axis-parallel boxes in dimension 1--6.
Intersection volumes and exact-membership atom volumes are computed by integer
Möbius inversion, without enumerating a coordinate grid.
"""
from __future__ import annotations

import itertools
import json
import math
import os
import random
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
DIMENSIONS = range(1, 7)
FAMILIES = ("common_core", "local_overlap", "r714_chain_random_axes", "r715_small_arms")
TOL = 1e-9
DENOMINATOR = 10**8
MIX = F(1, 10**12)


def enc(x):
    if not isinstance(x, F):
        x = F(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def box_vol(box):
    return math.prod(hi - lo for lo, hi in box)


def make_catalogue(d, family, local_index, seed):
    rng = random.Random(seed)
    count = 6 + (local_index % 2)
    boxes = []
    if family == "common_core":
        while len(boxes) < count:
            box = tuple((-rng.randint(1, 3), rng.randint(1, 4)) for _ in range(d))
            if box not in boxes:
                boxes.append(box)
    elif family == "local_overlap":
        for i in range(count):
            # The first-axis intervals form a chain of locally overlapping
            # neighbors but have empty total intersection.
            first = (2 * i, 2 * i + 3)
            rest = tuple((rng.randint(-2, 1), rng.randint(3, 6)) for _ in range(d - 1))
            box = (first,) + rest
            if box in boxes:
                # Change one other endpoint deterministically if needed.
                rest = tuple((lo, hi + (i + 1)) for lo, hi in rest)
                box = (first,) + rest
            boxes.append(box)
    elif family == "r714_chain_random_axes":
        xintervals = [(0, 10), (-2, 3), (7, 12), (2, 8), (5, 11),
                      (-1, 4), (6, 13)]
        for i in range(count):
            # A large interval, two separated end arms, and intermediate
            # chain members; all remaining coordinates randomly widen around
            # a fixed transverse core.
            first = xintervals[i]
            rest = tuple((-rng.randint(0, 2), 4 + rng.randint(0, 2)) for _ in range(d - 1))
            boxes.append((first,) + rest)
    elif family == "r715_small_arms":
        L = 4
        base = tuple((0, L) for _ in range(d))
        core = tuple((0, 1) if axis == 0 else (0, L) for axis in range(d))
        boxes = [base]
        ext = 1
        axis = 0
        while len(boxes) < count:
            row = list(core)
            lo, hi = row[axis]
            width = hi - lo
            row[axis] = (-ext * width, hi)
            candidate = tuple(row)
            if candidate not in boxes:
                boxes.append(candidate)
            ext += 1
            axis = (axis + 1) % d
    else:
        raise ValueError(f"unknown family {family}")

    assert len(boxes) == count and len(set(boxes)) == count
    assert all(len(box) == d and all(lo < hi and isinstance(lo, int) and isinstance(hi, int)
                                     for lo, hi in box) for box in boxes)
    return boxes


def intersection_volume(boxes, mask):
    indices = [i for i in range(len(boxes)) if mask & (1 << i)]
    if not indices:
        return 0
    volume = 1
    d = len(boxes[0])
    for axis in range(d):
        lo = max(boxes[i][axis][0] for i in indices)
        hi = min(boxes[i][axis][1] for i in indices)
        if hi <= lo:
            return 0
        volume *= hi - lo
    return volume


def exact_atoms(boxes):
    """Return positive exact-membership masks and their integer volumes."""
    nrows = len(boxes)
    maxmask = (1 << nrows) - 1
    intersections = {mask: intersection_volume(boxes, mask)
                     for mask in range(1, maxmask + 1)}
    exact = {}
    # Descending cardinality ensures all strict supersets have been computed.
    for mask in sorted(intersections, key=lambda v: (v.bit_count(), v), reverse=True):
        exact[mask] = intersections[mask] - sum(
            exact[sup] for sup in exact if sup != mask and (sup & mask) == mask)
        assert exact[mask] >= 0
    atoms = [(mask, volume) for mask, volume in exact.items() if volume > 0]
    union = sum(v for _, v in atoms)
    union_by_inclusion_exclusion = sum(
        ((-1) ** (mask.bit_count() + 1)) * vol
        for mask, vol in intersections.items())
    assert union == union_by_inclusion_exclusion
    # Independently verify each intersection as the sum of all containing atoms.
    for mask, vol in intersections.items():
        assert sum(v for atom_mask, v in atoms if (atom_mask & mask) == mask) == vol
    return atoms, union


def subfamily_atoms(atoms, selected_rows):
    merged = {}
    for oldmask, volume in atoms:
        newmask = 0
        for local_i, old_i in enumerate(selected_rows):
            if oldmask & (1 << old_i):
                newmask |= 1 << local_i
        if newmask:
            merged[newmask] = merged.get(newmask, 0) + volume
    result = sorted(merged.items())
    union = sum(v for _, v in result)
    # All selected rows must retain positive measure.
    for local_i in range(len(selected_rows)):
        bit = 1 << local_i
        assert sum(v for mask, v in result if mask & bit) > 0
    return result, union


def make_source_lp(box_volumes, atoms):
    r, M = len(box_volumes), len(atoms)
    A, b = [], []
    for i, volume in enumerate(box_volumes):
        A.append([F(mass) if mask & (1 << i) else F(0) for mask, mass in atoms])
        b.append(F(volume, 2))
    c = [F(mass) for _, mass in atoms]
    return A, b, c


def make_resp_lp(box_volumes, atoms, union):
    """Build the assigned lower-envelope LP stated in the r716 task."""
    r, M = len(box_volumes), len(atoms)
    A, b = [], []
    for j, (atommask, atomvol) in enumerate(atoms):
        for i in range(r):
            if atommask & (1 << i):
                row = [F(0)] * (2 * M)
                for k, (othermask, othervol) in enumerate(atoms):
                    if othermask & (1 << i):
                        row[k] = F(othervol, box_volumes[i])
                row[M + j] = F(-1)
                A.append(row)
                b.append(F(0))
    aggregate = [F(0)] * M + [F(mass) for _, mass in atoms]
    A.append(aggregate)
    b.append(F(union, 2))
    c = [F(mass) for _, mass in atoms] + [F(0)] * M
    # Strictly feasible reference used when repairing rationalized primals.
    reference = [F(1)] * M + [F(3, 4)] * M
    return A, b, c, reference


def _linprog(Aq, bq, cq, label):
    AA = np.asarray([[float(v) for v in row] for row in Aq], dtype=float)
    bb = np.asarray([float(v) for v in bq], dtype=float)
    cc = np.asarray([float(v) for v in cq], dtype=float)
    res = linprog(cc, A_ub=-AA, b_ub=-bb, bounds=(0, 1), method="highs",
                  options={"primal_feasibility_tolerance": TOL,
                           "dual_feasibility_tolerance": TOL})
    return res


def exact_certificate(A, b, c, result, label, reference=None, include_vectors=False):
    if reference is None:
        reference = [F(1)] * len(c)
    x0 = [min(F(1), max(F(0), F(float(v)).limit_denominator(DENOMINATOR)))
          for v in result.x]
    ref_slack = [sum((coef * rv for coef, rv in zip(row, reference)), F(0)) - rhs
                 for row, rhs in zip(A, b)]
    assert all(s > 0 for s in ref_slack), f"{label}: repair reference is not strictly feasible"
    eps = F(0)
    for row, rhs, rs in zip(A, b, ref_slack):
        slack = sum((coef * value for coef, value in zip(row, x0)), F(0)) - rhs
        if slack < 0:
            eps = max(eps, (-slack) / (rs - slack))
    eps += MIX
    assert eps < 1
    x = [(1 - eps) * v + eps * rv for v, rv in zip(x0, reference)]
    margins = [sum((coef * value for coef, value in zip(row, x)), F(0)) - rhs
               for row, rhs in zip(A, b)]
    assert all(s >= 0 for s in margins), f"{label}: repaired rational primal infeasible"
    assert all(F(0) <= v <= F(1) for v in x)
    primal = sum((ci * xi for ci, xi in zip(c, x)), F(0))

    y = [max(F(0), F(-float(v)).limit_denominator(DENOMINATOR))
         for v in result.ineqlin.marginals]
    coeff = [F(ci) for ci in c]
    for i, yi in enumerate(y):
        if yi:
            for j, aij in enumerate(A[i]):
                coeff[j] -= aij * yi
    dual = sum((bi * yi for bi, yi in zip(b, y)), F(0))
    dual += sum((min(F(0), v) for v in coeff), F(0))
    assert dual <= primal, f"{label}: exact weak duality failed"

    result_obj = {
        "solver_success": bool(result.success),
        "solver_status": int(result.status),
        "solver_message": str(result.message),
        "solver_objective_float": float(result.fun),
        "variables": len(c), "constraints": len(b),
        "tolerance": TOL,
        "rational_primal_upper": enc(primal),
        "rational_dual_lower": enc(dual),
        "objective_interval": [enc(dual), enc(primal)],
        "exact_relative_gap": enc((primal - dual) / dual) if dual > 0 else None,
        "primal_repair_mixture": enc(eps),
        "minimum_exact_constraint_slack": enc(min(margins)),
    }
    if include_vectors:
        result_obj["primal_rational"] = [enc(v) for v in x]
        result_obj["dual_nonzero_rational"] = [
            [i, enc(v)] for i, v in enumerate(y) if v]
    return result_obj


def certified_source_lp(box_volumes, atoms, include_vectors=False):
    A, b, c = make_source_lp(box_volumes, atoms)
    result = _linprog(A, b, c, "source")
    if not result.success:
        return {"solver_success": False, "solver_status": int(result.status),
                "solver_message": str(result.message), "tolerance": TOL}, None
    cert = exact_certificate(A, b, c, result, "source", include_vectors=include_vectors)
    return cert, result


def certified_resp_lp(box_volumes, atoms, union, include_vectors=False):
    A, b, c, ref = make_resp_lp(box_volumes, atoms, union)
    result = _linprog(A, b, c, "responsible")
    if not result.success:
        return {"solver_success": False, "solver_status": int(result.status),
                "solver_message": str(result.message), "tolerance": TOL}, None
    cert = exact_certificate(A, b, c, result, "responsible", reference=ref,
                             include_vectors=include_vectors)
    return cert, result


def scan_catalogue(task):
    d, family, local_index, seed = task
    boxes = make_catalogue(d, family, local_index, seed)
    original_atoms, original_union = exact_atoms(boxes)
    box_volumes = [box_vol(box) for box in boxes]
    candidate_rows = []
    failures = []
    best = None
    max_atoms_subfamily = 0

    for selected_mask in range(1, (1 << len(boxes))):
        selected = [i for i in range(len(boxes)) if selected_mask & (1 << i)]
        atoms, union = subfamily_atoms(original_atoms, selected)
        max_atoms_subfamily = max(max_atoms_subfamily, len(atoms))
        subvolumes = [box_volumes[i] for i in selected]
        source_cert, result = certified_source_lp(subvolumes, atoms, include_vectors=False)
        if result is None:
            failures.append({"selected_rows_0_based": selected,
                             "status": source_cert.get("solver_status"),
                             "message": source_cert.get("solver_message")})
            continue
        mlo, mhi = map(F, source_cert["objective_interval"])
        score_lo, score_hi = F(union, 1) / mhi, F(union, 1) / mlo
        score_est = union / float(result.fun)
        summary = {
            "selected_rows_0_based": selected,
            "selected_mask": selected_mask,
            "row_count": len(selected),
            "union_volume": union,
            "atom_count": len(atoms),
            "source_objective_interval": source_cert["objective_interval"],
            "score_interval": [enc(score_lo), enc(score_hi)],
            "score_estimate_float": score_est,
            "source_solver_status": source_cert["solver_status"],
        }
        candidate_rows.append(summary)
        if best is None or score_est > best["score_estimate_float"]:
            best = summary

    if best is None:
        return {"dimension": d, "family": family, "seed": seed,
                "boxes": boxes, "solver_failure": "all subfamily source LPs failed",
                "candidate_subfamily_count": (1 << len(boxes)) - 1,
                "solver_failures": failures}

    # Exact certificates for the numerically selected maximum-score subfamily.
    selected = best["selected_rows_0_based"]
    chosen_mask = best["selected_mask"]
    chosen_atoms, chosen_union = subfamily_atoms(original_atoms, selected)
    chosen_volumes = [box_volumes[i] for i in selected]
    source_cert, source_result = certified_source_lp(chosen_volumes, chosen_atoms,
                                                    include_vectors=True)
    resp_cert, resp_result = certified_resp_lp(chosen_volumes, chosen_atoms,
                                               chosen_union, include_vectors=True)
    if resp_result is None:
        resp_ratio = None
    else:
        mlo, mhi = map(F, source_cert["objective_interval"])
        rlo, rhi = map(F, resp_cert["objective_interval"])
        ratio_interval = [rlo / mhi, rhi / mlo]
        resp_ratio = {
            "numeric_estimate": float(resp_result.fun) / float(source_result.fun),
            "exact_certificate_interval": [enc(ratio_interval[0]), enc(ratio_interval[1])],
        }

    # Prove or report whether the selected candidate's score lower bound
    # separates it from every other subfamily's score upper bound.
    best_score_lower = F(best["score_interval"][0])
    competitor_upper = max((F(row["score_interval"][1]) for row in candidate_rows
                            if row["selected_mask"] != chosen_mask), default=F(0))
    certified_global_max = best_score_lower >= competitor_upper
    return {
        "dimension": d, "family": family, "seed": seed,
        "row_count": len(boxes), "boxes_integer_endpoints": boxes,
        "original_union_volume": original_union,
        "original_exact_atom_count": len(original_atoms),
        "candidate_subfamily_count": (1 << len(boxes)) - 1,
        "successful_source_lp_count": len(candidate_rows),
        "solver_failure_count": len(failures), "solver_failures": failures,
        "maximum_atoms_across_subfamilies": max_atoms_subfamily,
        "selected_max_score_subfamily": best,
        "global_max_separated_by_exact_intervals": certified_global_max,
        "max_other_subfamily_score_upper": enc(competitor_upper),
        "selected_score_lower": enc(best_score_lower),
        "source_LP_certificate": source_cert,
        "responsible_LP_certificate": resp_cert,
        "responsible_over_source_ratio": resp_ratio,
        "all_subfamily_screen": candidate_rows,
    }


def publication_summary(result):
    """Drop the exhaustive per-subset ledger and duplicate LP vectors."""
    selected = result["selected_max_score_subfamily"]
    return {
        "dimension": result["dimension"],
        "family": result["family"],
        "seed": result["seed"],
        "row_count": result["row_count"],
        "original_union_volume": result["original_union_volume"],
        "original_exact_atom_count": result["original_exact_atom_count"],
        "candidate_subfamily_count": result["candidate_subfamily_count"],
        "successful_source_lp_count": result["successful_source_lp_count"],
        "solver_failure_count": result["solver_failure_count"],
        "maximum_atoms_across_subfamilies": result["maximum_atoms_across_subfamilies"],
        "selected_max_score_subfamily": {
            "selected_rows_0_based": selected["selected_rows_0_based"],
            "selected_mask": selected["selected_mask"],
            "row_count": selected["row_count"],
            "union_volume": selected["union_volume"],
            "atom_count": selected["atom_count"],
            "source_objective_interval": selected["source_objective_interval"],
            "score_interval": selected["score_interval"],
            "score_estimate_float": selected["score_estimate_float"],
            "source_solver_status": selected["source_solver_status"],
        },
        "global_max_separated_by_exact_intervals": result["global_max_separated_by_exact_intervals"],
        "max_other_subfamily_score_upper": result["max_other_subfamily_score_upper"],
        "selected_score_lower": result["selected_score_lower"],
        "source_LP_status": {
            "success": result["source_LP_certificate"]["solver_success"],
            "status": result["source_LP_certificate"]["solver_status"],
            "tolerance": result["source_LP_certificate"]["tolerance"],
        },
        "responsible_LP_status": {
            "success": result["responsible_LP_certificate"]["solver_success"],
            "status": result["responsible_LP_certificate"]["solver_status"],
            "tolerance": result["responsible_LP_certificate"]["tolerance"],
        },
        "responsible_over_source_ratio": result["responsible_over_source_ratio"],
    }


def publication_view(full_output, full_archive):
    compact = dict(full_output)
    compact["full_detail_archive"] = full_archive
    compact["all_catalogues"] = [publication_summary(r) for r in full_output["all_catalogues"]]
    return compact


def main():
    tasks = []
    seed0 = 716202611
    case_i = 0
    for d in DIMENSIONS:
        for family_i, family in enumerate(FAMILIES):
            for local_index in range(5):
                seed = seed0 + case_i * 97 + d * 11 + family_i
                tasks.append((d, family, local_index, seed))
                case_i += 1

    start = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(scan_catalogue, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda r: (r["dimension"], FAMILIES.index(r["family"]), r["seed"]))
    valid = [r for r in results if "solver_failure" not in r and r["responsible_over_source_ratio"]]
    min5 = sorted(valid, key=lambda r: (
        r["responsible_over_source_ratio"]["numeric_estimate"], r["seed"]))[:5]
    min5_full = []
    for rank, result in enumerate(min5, 1):
        min5_full.append({
            "rank": rank,
            "dimension": result["dimension"], "family": result["family"],
            "seed": result["seed"], "boxes_integer_endpoints": result["boxes_integer_endpoints"],
            "selected_max_score_subfamily": result["selected_max_score_subfamily"],
            "global_max_separated_by_exact_intervals": result["global_max_separated_by_exact_intervals"],
            "source_LP_certificate": result["source_LP_certificate"],
            "responsible_LP_certificate": result["responsible_LP_certificate"],
            "responsible_over_source_ratio": result["responsible_over_source_ratio"],
        })

    output = {
        "status": "complete" if len(results) == len(tasks) else "incomplete",
        "scope_warning": "Finite random search diagnostic only; cannot prove a uniform constant or its failure.",
        "method": "Enumerate every nonempty row subset; exact integer intersection volumes and Möbius membership atoms; solve finite source LP per subset; solve responsible lower-envelope LP only for score-maximizing subset.",
        "catalogue_count": len(results), "target_catalogue_count": len(tasks),
        "dimensions": list(DIMENSIONS), "families": list(FAMILIES),
        "rows_per_catalogue": [6, 7],
        "subfamily_count_total": sum(r.get("candidate_subfamily_count", 0) for r in results),
        "successful_source_lp_count_total": sum(r.get("successful_source_lp_count", 0) for r in results),
        "solver_failure_count_total": sum(r.get("solver_failure_count", 0) for r in results),
        "catalogues_with_certified_global_score_max": sum(
            bool(r.get("global_max_separated_by_exact_intervals")) for r in results),
        "max_original_exact_atom_count": max(r.get("original_exact_atom_count", 0) for r in results),
        "max_subfamily_atom_count": max(r.get("maximum_atoms_across_subfamilies", 0) for r in results),
        "workers": 12, "blas_threads_per_worker": 1,
        "solver": "SciPy HiGHS linprog",
        "tolerance": {"primal_feasibility": TOL, "dual_feasibility": TOL},
        "rational_certificates": "Selected source and responsible LP primals are rationalized and repaired toward a strictly feasible rational reference; nonnegative rationalized dual multipliers yield exact weak-duality lower bounds.",
        "minimum_five_selected_candidates": min5_full,
        "all_catalogues": results,
        "wall_seconds": time.time() - start,
    }
    compact = publication_view(output,
        "/Volumes/data/CodexProjects/T/work/strong-maximal-drafts/r716/pruning_search_full.json")
    (HERE / "pruning_search.json").write_text(json.dumps(compact, indent=2) + "\n")
    print(json.dumps({k: output[k] for k in (
        "status", "catalogue_count", "subfamily_count_total",
        "successful_source_lp_count_total", "solver_failure_count_total",
        "catalogues_with_certified_global_score_max", "max_original_exact_atom_count",
        "max_subfamily_atom_count", "minimum_five_selected_candidates", "wall_seconds")}, indent=2))


if __name__ == "__main__":
    main()
