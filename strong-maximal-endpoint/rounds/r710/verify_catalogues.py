#!/usr/bin/env python3
"""Exact independent checker for the finite-catalogue certificates in catalogues.json.

No optimization is performed here. Geometry is reconstructed only by importing
catalogues.geometry; all certificate and objective arithmetic uses Fractions.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from catalogues import geometry  # permitted geometry constructor only


def q(value):
    """Parse an exact rational represented as an integer/fraction string."""
    if not isinstance(value, str):
        raise TypeError(f"expected rational string, got {value!r}")
    return Fraction(value)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def check_result(case, mode, result):
    label = f"n={case['n']} seed={case['seed']} mode={mode}"
    boxes = case["boxes"]
    common = bool(case["common_point"])
    Braw, massesraw, volraw, edges = geometry(boxes, common)
    B = [[bool(v) for v in row] for row in Braw.tolist()]
    masses = [int(v) for v in massesraw.tolist()]
    volumes = [int(v) for v in volraw.tolist()]
    M, r = len(masses), len(boxes)
    V = sum(masses)
    expected_edges = edges if mode == "isotone" else []
    cert = result["certificate"]

    require(result["atoms"] == M, f"{label}: atom count mismatch")
    require(result["rows"] == r, f"{label}: row count mismatch")
    require(result["union_volume"] == V, f"{label}: union volume mismatch")
    require(result["order_edges"] == len(expected_edges), f"{label}: edge count mismatch")
    require(all(sum(masses[j] for j in range(M) if B[i][j]) == volumes[i]
                for i in range(r)), f"{label}: box atom volumes do not match")

    # Recover the base primal vector. Omitted atom indices mean base value zero.
    base = [Fraction(0) for _ in range(M)]
    occupied = set()
    for value_s, indices in cert["base_primal_groups"].items():
        value = q(value_s)
        require(0 <= value <= 1, f"{label}: base primal value outside [0,1]")
        require(isinstance(indices, list), f"{label}: atom group is not a list")
        for j in indices:
            require(isinstance(j, int) and 0 <= j < M,
                    f"{label}: invalid base atom index {j!r}")
            require(j not in occupied, f"{label}: repeated base atom index {j}")
            occupied.add(j)
            base[j] = value

    eps = q(cert["strict_mix"])
    require(0 < eps < 1, f"{label}: strict mixture is not in (0,1)")
    if "strict_mixture" in result:
        require(q(result["strict_mixture"]) == eps,
                f"{label}: strict_mixture metadata mismatch")
    x = [(1 - eps) * value + eps for value in base]
    require(all(0 <= value <= 1 for value in x), f"{label}: mixed x outside [0,1]")
    primal = sum((masses[j] * x[j] for j in range(M)), Fraction(0))

    # Strict halo constraints, including all boxes (not a sampled subset).
    margins = []
    for i in range(r):
        box_mass = sum((masses[j] * x[j] for j in range(M) if B[i][j]), Fraction(0))
        margin = box_mass / volumes[i] - Fraction(1, 2)
        require(margin > 0, f"{label}: box {i} is not strictly above density 1/2")
        margins.append(margin)

    if mode == "isotone":
        require(common, f"{label}: isotone mode without common-point geometry")
        for k, (near, far) in enumerate(expected_edges):
            require(x[near] >= x[far],
                    f"{label}: isotonicity fails on edge {k} ({near},{far})")
    else:
        require(not cert["dual_order_multipliers"],
                f"{label}: plain certificate has order multipliers")

    # Dual multipliers are in source-mass units. Coefficients correspond to
    # c - B^T lambda - O^T eta with c_j = atom mass_j.
    lam = [Fraction(0) for _ in range(r)]
    seen_rows = set()
    for i, value_s in cert["dual_box_multipliers"]:
        require(isinstance(i, int) and 0 <= i < r,
                f"{label}: invalid dual box row {i!r}")
        require(i not in seen_rows, f"{label}: duplicate dual box row {i}")
        seen_rows.add(i)
        value = q(value_s)
        require(value >= 0, f"{label}: negative box multiplier")
        lam[i] = value

    eta = [Fraction(0) for _ in range(len(expected_edges))]
    seen_edges = set()
    for k, value_s in cert["dual_order_multipliers"]:
        require(mode == "isotone", f"{label}: order multiplier in plain mode")
        require(isinstance(k, int) and 0 <= k < len(expected_edges),
                f"{label}: invalid dual order edge {k!r}")
        require(k not in seen_edges, f"{label}: duplicate dual order edge {k}")
        seen_edges.add(k)
        value = q(value_s)
        require(value >= 0, f"{label}: negative order multiplier")
        eta[k] = value

    dual_coeff = [Fraction(m) for m in masses]
    for i, multiplier in enumerate(lam):
        if multiplier:
            for j in range(M):
                if B[i][j]:
                    dual_coeff[j] -= masses[j] * multiplier
    for multiplier, (near, far) in zip(eta, expected_edges):
        dual_coeff[near] -= multiplier
        dual_coeff[far] += multiplier
    dual = sum((lam[i] * volumes[i] / 2 for i in range(r)), Fraction(0))
    dual += sum((min(Fraction(0), c) for c in dual_coeff), Fraction(0))

    source_interval = [q(value) for value in result["source_interval"]]
    require(source_interval == [dual, primal],
            f"{label}: source interval mismatch; replay={[str(dual), str(primal)]}, "
            f"stored={[str(v) for v in source_interval]}")
    score_interval = [q(value) for value in result["fixed_catalogue_score_interval"]]
    expected_score = [Fraction(V, 1) / primal, Fraction(V, 1) / dual]
    require(score_interval == expected_score,
            f"{label}: fixed catalogue score interval mismatch")
    if "strict_source_mean_margin" in result:
        require(q(result["strict_source_mean_margin"]) == min(margins),
                f"{label}: strict source mean margin metadata mismatch")

    gap = (primal - dual) / dual
    require(gap >= 0, f"{label}: negative primal-dual gap")
    return {
        "n": case["n"],
        "K": case["K"],
        "seed": case["seed"],
        "common_point": common,
        "mode": mode,
        "primal": str(primal),
        "dual": str(dual),
        "relative_gap": str(gap),
        "min_box_margin": str(min(margins)),
        "atoms": M,
        "rows": r,
        "order_edges": len(expected_edges),
    }


def main():
    cases = json.loads((HERE / "catalogues.json").read_text())
    require(len(cases) == 40, f"expected 40 cases, found {len(cases)}")
    checked = []
    by_case = {}
    for case in cases:
        modes = [("plain", case["plain"])]
        if case["isotone"] is not None:
            modes.append(("isotone", case["isotone"]))
        for mode, result in modes:
            checked_result = check_result(case, mode, result)
            checked.append(checked_result)
            by_case[(case["n"], case["seed"])] = by_case.get((case["n"], case["seed"]), {})
            by_case[(case["n"], case["seed"])][mode] = result

    common_cases = [case for case in cases if case["common_point"]]
    intersections = 0
    for case in common_cases:
        plain_lo, plain_hi = map(q, case["plain"]["source_interval"])
        iso_lo, iso_hi = map(q, case["isotone"]["source_interval"])
        require(max(plain_lo, iso_lo) <= min(plain_hi, iso_hi),
                f"n={case['n']} seed={case['seed']}: plain/isotone intervals disjoint")
        intersections += 1

    gaps = [q(result["relative_gap"]) for result in checked]
    summary = {
        "status": "verified",
        "cases": len(cases),
        "certificates": len(checked),
        "plain_certificates": sum(1 for item in checked if item["mode"] == "plain"),
        "isotone_certificates": sum(1 for item in checked if item["mode"] == "isotone"),
        "common_point_interval_intersections": intersections,
        "all_primal_dual_certificate_gaps_nonnegative": all(g >= 0 for g in gaps),
        "maximum_exact_relative_gap": str(max(gaps, default=Fraction(0))),
        "maximum_relative_gap_certificate": checked[gaps.index(max(gaps))] if gaps else None,
        "certificates_detail": checked,
        "method": "exact Fraction arithmetic; geometry reconstructed via catalogues.geometry; no solver invoked",
    }
    (HERE / "replay.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: summary[key] for key in (
        "status", "cases", "certificates", "plain_certificates",
        "isotone_certificates", "common_point_interval_intersections",
        "maximum_exact_relative_gap", "maximum_relative_gap_certificate", "method")}, indent=2))


if __name__ == "__main__":
    main()
