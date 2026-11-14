#!/usr/bin/env python3
"""Priority-vertex enumeration and KL diagnostics for 12 rectangle catalogues.

All box intersections, exact-membership atoms, responsible weights, projection
prices, and TV right-hand sides use Fraction arithmetic.  KL values and local
strong-convexity left sides use mpmath at 100 decimal digits; they are numeric
diagnostics, not rigorous interval certificates.
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
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import mpmath as mp
import numpy as np
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
MP_DPS = 100
LP_TOL = 1e-9
QS = (F(1, 4), F(1, 2), F(3, 4))


def enc(x):
    if not isinstance(x, F):
        x = F(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def mpq(x):
    return mp.mpf(x.numerator) / x.denominator


def box_vol(box):
    return math.prod((hi - lo for lo, hi in box), start=F(1))


def intersection_volume(boxes, mask):
    ids = [i for i in range(len(boxes)) if mask & (1 << i)]
    if not ids:
        return F(0)
    d = len(boxes[0])
    vol = F(1)
    for axis in range(d):
        lo = max(boxes[i][axis][0] for i in ids)
        hi = min(boxes[i][axis][1] for i in ids)
        if hi <= lo:
            return F(0)
        vol *= hi - lo
    return vol


def membership_atoms(boxes):
    r = len(boxes)
    intersections = {mask: intersection_volume(boxes, mask)
                     for mask in range(1, 1 << r)}
    exact = {}
    for mask in sorted(intersections, key=lambda x: (x.bit_count(), x), reverse=True):
        exact[mask] = intersections[mask] - sum(
            value for sup, value in exact.items()
            if sup != mask and (sup & mask) == mask)
        assert exact[mask] >= 0
    atoms = sorted((mask, v) for mask, v in exact.items() if v > 0)
    for mask, value in intersections.items():
        assert sum((v for atom, v in atoms if (atom & mask) == mask), F(0)) == value
    union = sum((v for _, v in atoms), F(0))
    inc_exc = sum(((-1) ** (mask.bit_count() + 1)) * value
                  for mask, value in intersections.items())
    assert union == inc_exc
    return atoms, union


def build_inputs():
    base = ROOT
    wdata = json.loads((base / "rounds/r716/wolfram.json").read_text())["results"]
    unique = {}
    for record in wdata:
        key = json.dumps(record["boxes"], separators=(",", ":"))
        if key not in unique:
            unique[key] = {"boxes": record["boxes"], "dimension": record["dimension"]}
    assert len(unique) == 9
    r716 = []
    for i, item in enumerate(unique.values(), 1):
        boxes = [[[(F(str(lo)), F(str(hi))) for lo, hi in row] for row in item["boxes"]][j]
                 for j in range(len(item["boxes"]))]
        # Pull the published responsibility benchmark at q=1/2 for this box set.
        matching = next(x for x in wdata if x["boxes"] == item["boxes"] and x["p"] == "1/2")
        r716.append({"name": f"r716_catalogue_{i:02d}", "origin": "r716/wolfram.json",
                     "boxes": boxes, "known_mresp_half": F(matching["responsible_Hall_dual"]),
                     "known_source_half": F(matching["true_source"])})

    r714 = {"name": "r714_nearmax_M120_a11_b30", "origin": "r714/audit.json",
            "boxes": [
                [[F(0), F(120)]],
                [[F(-11), F(30)]],
                [[F(90), F(131)]],
            ], "known_mresp_half": F(2911, 52), "known_source_half": F(60)}
    more = []
    r715_data = json.loads((base / "rounds/r715/exact_certificates.json").read_text())
    for n in (3, 4):
        eps = F(6, n + 12)
        a = F(1, 2)
        boxes = []
        boxes.append([[(F(0), F(1)) for _ in range(n)]])
        for axis in range(n):
            row = [(F(0), F(1)) for _ in range(n)]
            if axis == 0:
                row[0] = (-a * eps, eps)
            else:
                row[0] = (F(0), eps)
                row[axis] = (-a, F(1))
            boxes.append([row])
        # Simplify nested singleton wrappers used above.
        boxes = [b[0] for b in boxes]
        # At a=1/2 and epsilon=epsilon*=6/(n+12), the r715 exact
        # certificate has b_{1/2}(P)=epsilon*.  Use that closed rational
        # value directly: r715's recorded dimensions need not include n=3.
        mresp = eps
        more.append({"name": f"r715_star_n{n}_a1over2", "origin": "r715/exact_certificates.json",
                     "boxes": boxes, "known_mresp_half": mresp,
                     "known_source_half": F(1, 2)})
    all_inputs = r716 + [r714] + more
    assert len(all_inputs) == 12
    for item in all_inputs:
        assert len(item["boxes"]) <= 6
    return all_inputs


def priority_vertices(atoms, row_vols, union):
    r = len(row_vols)
    pi_to_priority = {}
    permutation_count = 0
    for priority in itertools.permutations(range(r)):
        permutation_count += 1
        rank = [0] * r
        for position, row in enumerate(priority):
            rank[row] = position
        alloc = [F(0)] * r
        atom_owner = []
        for mask, vol in atoms:
            owner = min((i for i in range(r) if mask & (1 << i)), key=lambda i: rank[i])
            alloc[owner] += vol
            atom_owner.append(owner)
        pi = tuple(v / union for v in alloc)
        if pi not in pi_to_priority:
            pi_to_priority[pi] = {"priority": priority, "atom_owner": tuple(atom_owner)}
    assert all(sum(pi, F(0)) == 1 for pi in pi_to_priority)
    return pi_to_priority, permutation_count


def density_profile(boxes, atoms, union, pi):
    row_vols = [box_vol(box) for box in boxes]
    dens = []
    for mask, vol in atoms:
        d = sum((pi[i] / row_vols[i] for i in range(len(boxes)) if mask & (1 << i)), F(0))
        assert d > 0
        dens.append(d)
    assert sum((vol * d for (_, vol), d in zip(atoms, dens)), F(0)) == 1
    return dens


def kl_value(atoms, densities, union):
    value = mp.mpf("0")
    for (_, volume), density in zip(atoms, densities):
        mass = volume * density
        if mass:
            value += mpq(mass) * mp.log(mpq(density * union))
    return value


def np_price(atoms, densities, q):
    levels = sorted(set(densities), reverse=True)
    price = F(0)
    captured = F(0)
    for level in levels:
        group_volume = sum((vol for (_, vol), density in zip(atoms, densities)
                            if density == level), F(0))
        group_mass = group_volume * level
        if captured + group_mass < q:
            price += group_volume
            captured += group_mass
        else:
            price += (q - captured) / level
            captured = q
            break
    assert captured == q
    return price


def row_tv(boxes, i, j):
    inter = intersection_volume(boxes, (1 << i) | (1 << j))
    vi, vj = box_vol(boxes[i]), box_vol(boxes[j])
    return 1 - inter / max(vi, vj)


def local_switch_inequality(boxes, atoms, union, densities, atom_owner, pi):
    row_vols = [box_vol(b) for b in boxes]
    ell = []
    for i, vol_i in enumerate(row_vols):
        value = mp.mpf("0")
        for (mask, atomvol), density in zip(atoms, densities):
            if mask & (1 << i):
                value += mpq(atomvol / vol_i) * mp.log(mpq(density))
        ell.append(value)
    worst_margin = None
    failed = []
    checked = 0
    for (mask, atomvol), owner in zip(atoms, atom_owner):
        qatom = atomvol / union
        for j in range(len(boxes)):
            if j == owner or not (mask & (1 << j)):
                continue
            rhs = 2 * qatom * row_tv(boxes, owner, j)**2
            lhs = ell[owner] - ell[j]
            margin = lhs - mpq(rhs)
            checked += 1
            if worst_margin is None or margin < worst_margin:
                worst_margin = margin
                tight = {"owner": owner, "alternative": j,
                         "atom_membership_mask": mask,
                         "atom_volume": enc(atomvol), "Q_atom": enc(qatom),
                         "TV_exact": enc(row_tv(boxes, owner, j)),
                         "rhs_exact": enc(rhs), "lhs_numeric": mp.nstr(lhs, MP_DPS),
                         "margin_numeric": mp.nstr(margin, MP_DPS)}
            if margin < mp.mpf("-1e-85"):
                failed.append({"owner": owner, "alternative": j,
                               "atom_mask": mask,
                               "margin_numeric": mp.nstr(margin, MP_DPS)})
    return {"checked_inequalities": checked, "all_pass_numerically": not failed,
            "numeric_tolerance": "1e-85", "worst_margin": tight if checked else None,
            "failure_count": len(failed), "first_failures": failed[:5]}


def globalmaxscore_diagnostic(boxes):
    """Float LP scan over all subfamilies; signed margins are diagnostic only."""
    r = len(boxes)
    records = []
    for subset in range(1, 1 << r):
        selected = [i for i in range(r) if subset & (1 << i)]
        subboxes = [boxes[i] for i in selected]
        atoms, union = membership_atoms(subboxes)
        vols = [box_vol(b) for b in subboxes]
        A = [[float(v if (mask & (1 << i)) else 0)
              for mask, v in atoms] for i in range(len(subboxes))]
        b = [-float(v) / 2 for v in vols]
        c = [float(v) for _, v in atoms]
        result = linprog(c, A_ub=-np.asarray(A), b_ub=np.asarray(b),
                         bounds=(0, 1), method="highs",
                         options={"primal_feasibility_tolerance": LP_TOL,
                                  "dual_feasibility_tolerance": LP_TOL})
        if not result.success:
            continue
        score = float(union) / float(result.fun)
        records.append({"subset": selected, "rows": len(selected),
                        "union": float(union), "source_mass_float": float(result.fun),
                        "score_float": score, "status": int(result.status)})
    full_mask = (1 << r) - 1
    full = next(x for x in records if len(x["subset"]) == r)
    proper_best = max((x for x in records if len(x["subset"]) < r),
                      key=lambda x: x["score_float"], default=None)
    global_best = max(records, key=lambda x: x["score_float"])
    margin_full_vs_proper = (full["score_float"] - proper_best["score_float"]
                             if proper_best else None)
    return {
        "subfamilies_checked": len(records), "solver_failures": (1 << r) - 1 - len(records),
        "full_catalogue_score_float": full["score_float"],
        "best_proper_subcatalogue": proper_best,
        "globalmaxscore_float": global_best["score_float"],
        "globalmaxscore_rows": global_best["subset"],
        "globalmaxscore_is_full_catalogue": len(global_best["subset"]) == r,
        "full_minus_best_proper_signed_margin_float": margin_full_vs_proper,
        "full_minus_globalmax_signed_margin_float": full["score_float"] - global_best["score_float"],
        "sign_status": ("positive" if margin_full_vs_proper > 0 else
                        "negative" if margin_full_vs_proper < 0 else "zero") if proper_best else "no_proper_subfamily",
        "solver_method": "SciPy HiGHS float LP, tolerance 1e-9; numerical diagnostic only, not 100-digit strict interval certification",
    }


def run_case(case):
    mp.mp.dps = MP_DPS
    boxes = case["boxes"]
    atoms, union = membership_atoms(boxes)
    row_vols = [box_vol(b) for b in boxes]
    vertices, permutations = priority_vertices(atoms, row_vols, union)
    scored = []
    per_q = {q: [] for q in QS}
    for pi, witness in vertices.items():
        dens = density_profile(boxes, atoms, union, pi)
        kl = kl_value(atoms, dens, union)
        prices = {q: np_price(atoms, dens, q) for q in QS}
        for q in QS:
            per_q[q].append((prices[q], pi, witness))
        scored.append((kl, pi, witness, dens, prices))
    best = max(scored, key=lambda x: x[0])
    kl, best_pi, witness, densities, best_prices = best
    prices_max = {}
    for q in QS:
        max_price = max(item[0] for item in per_q[q])
        price_witness = next(item for item in per_q[q] if item[0] == max_price)
        prices_max[str(float(q))] = {
            "q": enc(q), "max_bq_over_priority_vertices": enc(max_price),
            "maximizing_pi": [enc(x) for x in price_witness[1]],
        }

    comparison = None
    if case.get("known_mresp_half") is not None:
        computed = F(prices_max["0.5"]["max_bq_over_priority_vertices"])
        known = case["known_mresp_half"]
        comparison = {"known_r716_or_formula": enc(known),
                     "computed_max_b_half": enc(computed),
                     "exact_match": computed == known}

    best_price_table = {str(float(q)): enc(best_prices[q]) for q in QS}
    switch = local_switch_inequality(boxes, atoms, union, densities,
                                     witness["atom_owner"], best_pi)
    global_diag = globalmaxscore_diagnostic(boxes)
    runner_up = max((x[0] for x in scored if x[1] != best_pi), default=kl)
    return {
        "name": case["name"], "origin": case["origin"],
        "dimension": len(boxes[0]), "row_count": len(boxes),
        "boxes_rational_endpoints": [[[enc(lo), enc(hi)] for lo, hi in box]
                                      for box in boxes],
        "union_volume": enc(union), "row_volumes": [enc(v) for v in row_vols],
        "membership_atom_count": len(atoms),
        "priority_permutations_enumerated": permutations,
        "unique_responsible_vertices": len(vertices),
        "priority_vertex_completeness_basis": "Coverage polymatroid base vertices are generated by all row-priority greedy orders; vectors deduplicated exactly.",
        "max_KL_vertex": {
            "priority_1_based_high_to_low": [i + 1 for i in witness["priority"]],
            "pi_rational": [enc(v) for v in best_pi],
            "atom_owner_1_based": [i + 1 for i in witness["atom_owner"]],
            "atom_density_rational": [enc(v) for v in densities],
            "KL_P_vs_uniform_Q_numeric_100dps": mp.nstr(kl, MP_DPS),
            "KL_runner_up_gap_numeric_100dps": mp.nstr(kl - runner_up, MP_DPS),
            "prices_bq_for_max_KL_pi": best_price_table,
        },
        "max_NP_price_over_all_responsible_vertices": prices_max,
        "known_responsible_half_price_match": comparison,
        "local_single_atom_switch_inequality_at_max_KL_vertex": switch,
        "globalmaxscore_diagnostic": global_diag,
    }


def main():
    mp.mp.dps = MP_DPS
    cases = build_inputs()
    start = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(run_case, case) for case in cases]
        for future in as_completed(futures):
            results.append(future.result())
    by_name = {r["name"]: r for r in results}
    results = [by_name[c["name"]] for c in cases]
    output = {
        "status": "complete" if len(results) == 12 else "incomplete",
        "purpose": "Finite diagnostic for max-KL responsible vertices, exact NP prices, and one-atom switch inequalities; not a general theorem.",
        "inputs": "Nine distinct r716/wolfram.json box sets, r714 nearmax (M=120,a=11,b=30), and r715 n=3,4 star-arm boxes with a=1/2, epsilon=6/(n+12).",
        "vertex_method": "Enumerated every priority permutation; each exact atom assigned to its first highest-priority covering row; duplicate rational pi vectors removed.",
        "kl_precision_decimal_digits": MP_DPS,
        "numerical_certification_note": "KL and switch-inequality left sides are 100-digit mpmath evaluations, not rigorous intervals. Global-max-score margins are float HiGHS diagnostics. Exact rational TV, atom volumes, pi, NP prices, and right sides are replayable.",
        "processes": 12, "blas_threads_per_worker": 1,
        "globalmaxscore_sign_summary": {
            "positive_full_minus_best_proper": sum(1 for x in results if x["globalmaxscore_diagnostic"]["sign_status"] == "positive"),
            "negative_full_minus_best_proper": sum(1 for x in results if x["globalmaxscore_diagnostic"]["sign_status"] == "negative"),
            "zero_full_minus_best_proper": sum(1 for x in results if x["globalmaxscore_diagnostic"]["sign_status"] == "zero"),
        },
        "known_half_price_matches": sum(1 for x in results
                                         if x["known_responsible_half_price_match"] and
                                         x["known_responsible_half_price_match"]["exact_match"]),
        "case_count": len(results), "results": results,
        "wall_seconds": time.time() - start,
    }
    (HERE / "greedy.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: output[k] for k in (
        "status", "case_count", "known_half_price_matches",
        "globalmaxscore_sign_summary", "wall_seconds")}, indent=2))


if __name__ == "__main__":
    main()
