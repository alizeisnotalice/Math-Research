#!/usr/bin/env python3
"""Local responsible-partition search for common-core star products.

For each capped (h,k) pair, sample 200 genuine atom-to-row receiving
partitions.  All geometry, responsible weights, atom densities, and projection
cap checks use exact Fractions.  Entropies use mpmath at 80 decimal digits and
are numerical diagnostics, not rigorous numerical intervals.
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

import mpmath as mp

HERE = Path(__file__).resolve().parent
MP_DPS = 80
PARTITIONS_PER_PARAMETER = 200
SEED0 = 719202611
CHOSEN_PAIRS = [(h, k) for h in (1, 2, 3, 4) for k in (2, 3, 4)
                if k**h <= 64 and (k + 1)**h <= 625]


def enc(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def mpq(x):
    return mp.mpf(x.numerator) / x.denominator


def make_structure(h, k):
    a = F(1, 2 * h)
    rows = list(itertools.product(range(k), repeat=h))
    atoms = []
    for mask in range(1 << h):
        coords = tuple(j for j in range(h) if mask & (1 << j))
        for partial in itertools.product(range(k), repeat=len(coords)):
            labels = dict(zip(coords, partial))
            covers = [i for i, row in enumerate(rows)
                      if all(row[j] == labels[j] for j in coords)]
            assert covers
            volume = a ** len(coords)
            atoms.append({"mask": mask, "coords": coords, "partial": partial,
                          "volume": volume, "covers": covers})
    omega = (1 + k * a) ** h
    v = (1 + a) ** h
    assert len(rows) == k**h and len(atoms) == (k + 1)**h
    assert sum((atom["volume"] for atom in atoms), F(0)) == omega
    assert v <= 2
    denom = (2 * h)**h
    atom_units = [int(atom["volume"] * denom) for atom in atoms]
    omega_units = int(omega * denom)
    assert all(x > 0 for x in atom_units) and sum(atom_units) == omega_units
    return a, rows, atoms, v, omega, atom_units, omega_units


def receiving_assignment(atoms, rows, atom_units, rng, method):
    assigned = []
    row_units = [0] * len(rows)
    atom_order = list(range(len(atoms)))
    if method == "greedy_concentrate":
        atom_order.sort(key=lambda i: (-atom_units[i], i))
    for atom_i in atom_order:
        atom = atoms[atom_i]
        covers = atom["covers"]
        if method == "canonical":
            labels = dict(zip(atom["coords"], atom["partial"]))
            row = tuple(labels.get(j, 0) for j in range(len(rows[0])))
            recipient = rows.index(row)
            assert recipient in covers
        elif method == "greedy_concentrate":
            recipient = max(covers, key=lambda i: (row_units[i], -i))
        else:
            recipient = rng.choice(covers)
        assigned.append((atom_i, recipient))
        row_units[recipient] += atom_units[atom_i]
    assert sorted(i for i, _ in assigned) == list(range(len(atoms)))
    assert sum(row_units) == sum(atom_units)
    return assigned, row_units


def entropy_of_probs(probs):
    result = mp.mpf("0")
    for probability in probs:
        if probability:
            p = mpq(probability)
            result -= p * mp.log(p)
    return result


def source_density_and_entropy(a, rows, atoms, pi, v):
    densities = []
    h_source = mp.mpf("0")
    for atom in atoms:
        mass = sum((pi[i] for i in atom["covers"]), F(0))
        density = mass / v
        densities.append(density)
        assert density >= 0
        if density:
            h_source -= mpq(atom["volume"] * density) * mp.log(mpq(density))
    total = sum((atom["volume"] * density
                 for atom, density in zip(atoms, densities)), F(0))
    assert total == 1
    return densities, h_source


def expected_projection_entropy(h, k, a, rows, pi):
    pext = a / (1 + a)
    expectation = mp.mpf("0")
    marginal_cache = {}
    for mask in range(1 << h):
        coords = tuple(j for j in range(h) if mask & (1 << j))
        counts = {}
        for row, weight in zip(rows, pi):
            label = tuple(row[j] for j in coords)
            counts[label] = counts.get(label, F(0)) + weight
        H = entropy_of_probs(counts.values())
        probability_J = mpq(pext ** len(coords) * (1 - pext) ** (h - len(coords)))
        expectation += probability_J * H
        marginal_cache[mask] = counts
    return expectation, marginal_cache


def verify_projection_caps(h, k, a, rows, pi):
    cap_base = (1 + a) / (1 + k * a)
    constraints = 0
    min_slack = None
    tightest = None
    for mask in range(1, 1 << h):
        coords = tuple(j for j in range(h) if mask & (1 << j))
        cap = cap_base ** len(coords)
        for partial in itertools.product(range(k), repeat=len(coords)):
            actual = sum((weight for row, weight in zip(rows, pi)
                          if all(row[j] == val for j, val in zip(coords, partial))), F(0))
            slack = cap - actual
            assert slack >= 0
            constraints += 1
            if min_slack is None or slack < min_slack:
                min_slack = slack
                tightest = {"coordinates_1_based": [j + 1 for j in coords],
                            "labels_1_based": [t + 1 for t in partial],
                            "actual": enc(actual), "cap": enc(cap), "slack": enc(slack)}
    assert constraints == (k + 1)**h - 1
    return {"verified": True, "constraint_count": constraints,
            "minimum_exact_slack": enc(min_slack), "tightest_constraint": tightest}


def evaluate_partition(h, k, a, rows, atoms, v, omega, atom_units,
                       omega_units, assigned, row_units, index, method, seed):
    pi = [F(units, omega_units) for units in row_units]
    assert sum(pi, F(0)) == 1
    densities, h_source = source_density_and_entropy(a, rows, atoms, pi, v)
    expected_H, _ = expected_projection_entropy(h, k, a, rows, pi)
    gap_direct = h_source - mp.log(mpq(v / 2))
    gap_formula = mp.log(2) + expected_H
    identity_abs_error = abs(gap_direct - gap_formula)
    assert identity_abs_error < mp.mpf("1e-60")

    pext = a / (1 + a)
    log_lower = mp.log(2) + h * mpq(pext) * mp.log(mpq((1 + k * a) / (1 + a)))
    log_upper = mp.log(2) + h * mpq(pext) * mp.log(k)
    margin_tol = mp.mpf("1e-60")
    assert gap_formula + margin_tol >= log_lower
    assert gap_formula <= log_upper + margin_tol
    cap_check = verify_projection_caps(h, k, a, rows, pi)

    return {
        "partition_index": index,
        "partition_seed": seed,
        "method": method,
        "gap_numeric": mp.nstr(gap_formula, MP_DPS),
        "h_source_numeric": mp.nstr(h_source, MP_DPS),
        "expected_projection_entropy_numeric": mp.nstr(expected_H, MP_DPS),
        "identity_absolute_error_numeric": mp.nstr(identity_abs_error, 8),
        "lower_bound_numeric": mp.nstr(log_lower, MP_DPS),
        "upper_bound_numeric": mp.nstr(log_upper, MP_DPS),
        "gap_within_bound": True,
        "projection_caps": cap_check,
        "pi": pi,
        "densities": densities,
        "assigned_rows": [recipient for _, recipient in sorted(assigned)],
    }


def run_parameter(job):
    h, k, parameter_index = job
    mp.mp.dps = MP_DPS
    a, rows, atoms, v, omega, atom_units, omega_units = make_structure(h, k)
    seed_base = SEED0 + parameter_index * 100003 + h * 1009 + k * 101
    minimum = None
    gap_values = []
    for index in range(PARTITIONS_PER_PARAMETER):
        seed = seed_base + index
        if index == 0:
            method = "canonical_default_completion"
            rng = random.Random(seed)
        elif index == 1:
            method = "greedy_concentrate"
            rng = random.Random(seed)
        else:
            method = "independent_uniform_covering_row_per_atom"
            rng = random.Random(seed)
        assigned, row_units = receiving_assignment(atoms, rows, atom_units, rng, method)
        evaluated = evaluate_partition(h, k, a, rows, atoms, v, omega,
                                       atom_units, omega_units, assigned,
                                       row_units, index, method, seed)
        gap = mp.mpf(evaluated["gap_numeric"])
        gap_values.append(gap)
        if minimum is None or gap < minimum[0]:
            minimum = (gap, evaluated)

    best = minimum[1]
    return {
        "h": h, "k": k, "row_count": len(rows), "atom_count": len(atoms),
        "a": enc(a), "v_row_volume": enc(v), "Omega_union_volume": enc(omega),
        "v_le_2": v <= 2,
        "seed_base": seed_base,
        "partition_seed_rule": "seed_base + partition_index",
        "partition_algorithms": {
            "index_0": "complete every unspecified row label with 1 (canonical)",
            "index_1": "greedy concentrate mass on an already-heavy covering row",
            "indices_2_to_199": "independent uniform choice among covering rows per atom using random.Random(seed)",
        },
        "sample_count": PARTITIONS_PER_PARAMETER,
        "all_sample_gaps_within_bounds": True,
        "all_sample_projection_caps_exact": True,
        "sampled_gap_min_numeric": mp.nstr(min(gap_values), MP_DPS),
        "sampled_gap_max_numeric": mp.nstr(max(gap_values), MP_DPS),
        "minimum_candidate": {
            "partition_index": best["partition_index"],
            "partition_seed": best["partition_seed"],
            "method": best["method"],
            "gap_numeric": best["gap_numeric"],
            "h_source_numeric": best["h_source_numeric"],
            "expected_projection_entropy_numeric": best["expected_projection_entropy_numeric"],
            "identity_absolute_error_numeric": best["identity_absolute_error_numeric"],
            "lower_bound_numeric": best["lower_bound_numeric"],
            "upper_bound_numeric": best["upper_bound_numeric"],
            "pi_rational": [enc(p) for p in best["pi"]],
            "density_per_atom_rational": [enc(d) for d in best["densities"]],
            "assigned_row_1_based_per_atom": [i + 1 for i in best["assigned_rows"]],
            "projection_caps_exact": best["projection_caps"],
        },
    }


def main():
    tasks = [(h, k, idx) for idx, (h, k) in enumerate(CHOSEN_PAIRS)]
    start = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(run_parameter, job) for job in tasks]
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda r: (r["h"], r["k"]))
    output = {
        "status": "complete",
        "scope": "Finite local star-product partition search only; the minimum among 200 sampled partitions is not a global optimum over responsible transports.",
        "source_draft": "/Volumes/data/CodexProjects/T/work/strong-maximal-drafts/r719/testing_entropy.md",
        "h_values_requested": [1, 2, 3, 4], "k_values_requested": [2, 3, 4],
        "excluded_by_row_cap": [{"h": h, "k": k, "row_count": k**h,
                                  "reason": "k^h exceeds the stated 64-row cap"}
                                 for h, k in itertools.product((1, 2, 3, 4), (2, 3, 4))
                                 if (h, k) not in CHOSEN_PAIRS],
        "parameters_run": len(results),
        "partitions_per_parameter": PARTITIONS_PER_PARAMETER,
        "partition_count_total": sum(r["sample_count"] for r in results),
        "workers": 12,
        "blas_threads_per_worker": 1,
        "mpmath_decimal_precision": MP_DPS,
        "gap_identity": "h(P)-log(v/2)=log(2)+E_J H(T_J)",
        "analytic_lower_bound": "log(2)+h*a/(1+a)*log((1+k*a)/(1+a))",
        "analytic_upper_bound": "log(2)+h*a/(1+a)*log(k)",
        "fraction_checks": "exact geometry/partition pi/cell densities/projection caps; entropy uses mpmath 80 dps and identity/bounds tolerance 1e-60",
        "all_parameters_passed": all(r["all_sample_gaps_within_bounds"] and
                                      r["all_sample_projection_caps_exact"] for r in results),
        "results": results,
        "wall_seconds": time.time() - start,
    }
    (HERE / "projection_search.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: output[k] for k in (
        "status", "parameters_run", "partitions_per_parameter",
        "partition_count_total", "workers", "blas_threads_per_worker",
        "mpmath_decimal_precision", "all_parameters_passed", "wall_seconds")}, indent=2))


if __name__ == "__main__":
    main()
