#!/usr/bin/env python3
"""Independent Fraction-only replay of the five published r716 certificates.

No solver is imported or called.  Atom membership is rebuilt from integer box
intersection volumes using exact Möbius inversion, then every primal/dual LP
constraint and objective bound is checked afresh.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent


def Q(x):
    return F(x)


def intersection_volume(boxes, mask):
    ids = [i for i in range(len(boxes)) if mask & (1 << i)]
    d = len(boxes[0])
    vol = 1
    for axis in range(d):
        lo = max(boxes[i][axis][0] for i in ids)
        hi = min(boxes[i][axis][1] for i in ids)
        if hi <= lo:
            return 0
        vol *= hi - lo
    return vol


def mobius_atoms(boxes):
    r = len(boxes)
    intersections = {m: intersection_volume(boxes, m) for m in range(1, 1 << r)}
    exact = {}
    for mask in sorted(intersections, key=lambda m: (m.bit_count(), m), reverse=True):
        exact[mask] = intersections[mask] - sum(
            v for sup, v in exact.items() if (sup & mask) == mask and sup != mask)
        assert exact[mask] >= 0
    atoms = sorted((mask, v) for mask, v in exact.items() if v > 0)
    for mask, value in intersections.items():
        assert sum(v for atom, v in atoms if (atom & mask) == mask) == value
    union = sum(v for _, v in atoms)
    ie = sum(((-1) ** (mask.bit_count() + 1)) * value
             for mask, value in intersections.items())
    assert union == ie
    return atoms, union


def dense_dual(cert, nrows):
    y = [F(0)] * nrows
    for index, value in cert["dual_nonzero_rational"]:
        assert 0 <= index < nrows and y[index] == 0
        y[index] = Q(value)
    assert all(v >= 0 for v in y)
    return y


def replay_lp(A, b, c, cert, label):
    nvar = len(c)
    assert cert["solver_success"] and cert["solver_status"] == 0
    assert cert["variables"] == nvar and cert["constraints"] == len(b)
    x = list(map(Q, cert["primal_rational"]))
    assert len(x) == nvar and all(0 <= v <= 1 for v in x), f"{label}: primal bounds"
    margins = [sum((a * v for a, v in zip(row, x)), F(0)) - rhs
               for row, rhs in zip(A, b)]
    assert all(v >= 0 for v in margins), f"{label}: primal infeasible"
    primal = sum((ci * xi for ci, xi in zip(c, x)), F(0))
    assert primal == Q(cert["rational_primal_upper"]), f"{label}: primal objective mismatch"

    y = dense_dual(cert, len(b))
    coeff = list(map(F, c))
    for i, yi in enumerate(y):
        for j, aij in enumerate(A[i]):
            coeff[j] -= aij * yi
    dual = sum((bi * yi for bi, yi in zip(b, y)), F(0))
    dual += sum((min(F(0), v) for v in coeff), F(0))
    assert dual == Q(cert["rational_dual_lower"]), f"{label}: dual objective mismatch"
    assert dual <= primal, f"{label}: weak duality"
    assert [str(dual) if dual.denominator != 1 else str(dual.numerator),
            str(primal) if primal.denominator != 1 else str(primal.numerator)] == cert["objective_interval"]
    assert min(margins) == Q(cert["minimum_exact_constraint_slack"])
    assert (primal - dual) / dual == Q(cert["exact_relative_gap"])
    return {"primal_upper": str(primal), "dual_lower": str(dual),
            "min_constraint_slack": str(min(margins)),
            "relative_gap": str((primal - dual) / dual),
            "variables": nvar, "constraints": len(b)}


def source_lp(boxes, atoms):
    vols = [math.prod(hi - lo for lo, hi in box) for box in boxes]
    A = [[F(mass) if mask & (1 << i) else F(0) for mask, mass in atoms]
         for i in range(len(boxes))]
    b = [F(v, 2) for v in vols]
    c = [F(mass) for _, mass in atoms]
    return A, b, c


def responsible_lp(boxes, atoms, union):
    r, M = len(boxes), len(atoms)
    vols = [math.prod(hi - lo for lo, hi in box) for box in boxes]
    A, b = [], []
    for j, (atommask, atomvol) in enumerate(atoms):
        for i in range(r):
            if atommask & (1 << i):
                row = [F(0)] * (2 * M)
                for k, (othermask, othervol) in enumerate(atoms):
                    if othermask & (1 << i):
                        row[k] = F(othervol, vols[i])
                row[M + j] = F(-1)
                A.append(row)
                b.append(F(0))
    A.append([F(0)] * M + [F(v) for _, v in atoms])
    b.append(F(union, 2))
    c = [F(v) for _, v in atoms] + [F(0)] * M
    return A, b, c


def encode(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def replay_candidate(candidate):
    selected = candidate["selected_max_score_subfamily"]["selected_rows_0_based"]
    all_boxes = candidate["boxes_integer_endpoints"]
    boxes = [[tuple(map(int, axis)) for axis in all_boxes[i]] for i in selected]
    atoms, union = mobius_atoms(boxes)
    meta = candidate["selected_max_score_subfamily"]
    assert union == meta["union_volume"]
    assert len(atoms) == meta["atom_count"]

    A_s, b_s, c_s = source_lp(boxes, atoms)
    source_result = replay_lp(A_s, b_s, c_s, candidate["source_LP_certificate"],
                              "source")
    A_r, b_r, c_r = responsible_lp(boxes, atoms, union)
    resp_result = replay_lp(A_r, b_r, c_r, candidate["responsible_LP_certificate"],
                            "responsible")

    source_interval = [Q(v) for v in candidate["source_LP_certificate"]["objective_interval"]]
    resp_interval = [Q(v) for v in candidate["responsible_LP_certificate"]["objective_interval"]]
    ratio_interval = [resp_interval[0] / source_interval[1],
                      resp_interval[1] / source_interval[0]]
    assert [encode(v) for v in ratio_interval] == candidate["responsible_over_source_ratio"]["exact_certificate_interval"]
    score_interval = [F(union) / source_interval[1], F(union) / source_interval[0]]
    assert [encode(v) for v in score_interval] == meta["score_interval"]
    return {
        "rank": candidate["rank"], "dimension": candidate["dimension"],
        "family": candidate["family"], "seed": candidate["seed"],
        "selected_rows_0_based": selected,
        "union_volume": union, "mobius_atom_count": len(atoms),
        "all_intersection_volumes_reconstructed": True,
        "source_LP": source_result, "responsible_LP": resp_result,
        "ratio_certificate_interval": [encode(v) for v in ratio_interval],
        "score_certificate_interval": [encode(v) for v in score_interval],
        "status": "verified",
    }


def main():
    published = json.loads((HERE / "pruning_search.json").read_text())
    assert len(published["minimum_five_selected_candidates"]) == 5
    replayed = [replay_candidate(c) for c in published["minimum_five_selected_candidates"]]
    out = {
        "status": "all_five_exact_replays_passed",
        "method": "Independent Python Fraction reconstruction; no optimizer or solver imported",
        "candidate_count": len(replayed),
        "mobius_intersection_reconstruction_passes": len(replayed),
        "source_primal_dual_replay_passes": len(replayed),
        "responsible_primal_dual_replay_passes": len(replayed),
        "results": replayed,
    }
    (HERE / "replay.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in (
        "status", "candidate_count", "mobius_intersection_reconstruction_passes",
        "source_primal_dual_replay_passes", "responsible_primal_dual_replay_passes")}, indent=2))


if __name__ == "__main__":
    main()
