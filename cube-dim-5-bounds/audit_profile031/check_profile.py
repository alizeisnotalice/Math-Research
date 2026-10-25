#!/usr/bin/env python3
"""Finite independent scalar/eigenspace audit of the continuous-rho patch barrier."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plus(a, b):
    return [Q(a[k] if k < len(a) else 0)+Q(b[k] if k < len(b) else 0)
            for k in range(max(len(a), len(b)))]


def times(a, b):
    result = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def scale(a, x):
    return [Q(x)*y for y in a]


def evaluate(sources):
    saved = json.loads(sources["candidate_result"].read_text())
    # Four head-sign patterns only, followed by the four free tail choices.
    head_counts = Counter((1+a+b)**2 for a, b in product((-1, 1), repeat=2))
    direction_counts = {str(k): 4*v for k, v in head_counts.items()}
    assert direction_counts == saved["direction_counts"] == {"9": 4, "1": 12}
    # Eigenspaces: head parallel, two head-orthogonal directions, two tails.
    parallel, head_orthogonal, tail = [4, 12], [16], [12, 4]
    determinant = times(times(times(times(parallel, head_orthogonal), head_orthogonal), tail), tail)
    bD = scale(times(times(head_orthogonal, head_orthogonal), times(tail, tail)), 3)
    assert determinant == scale(times([1, 3], times([3, 1], [3, 1])), 16384)
    assert bD == scale(times([3, 1], [3, 1]), 12288)
    # Rational identities for all rho<1, rather than any sampled rho values.
    numerator, denominator = [12, -1], [12, -5]
    lower_difference = plus(scale(numerator, 5), scale(denominator, -1))
    upper_difference = plus(scale(denominator, 11), scale(numerator, -7))
    assert lower_difference == [48, 0]
    assert upper_difference == [48, -48]
    assert Q(11, 7) < Q(9, 7)**2 and Q(9, 7) < Q(8, 7)**2
    derivative_numerator = plus(scale([3, 1], 3), scale([1, 3], -1))
    assert derivative_numerator == [8, 0]
    assert Q(2, 3) < Q(5, 6)**2
    # Unshifted positive Mellin integrals after y=C sqrt(z).
    gamma9half_over_sqrtpi = Q(1)
    for k in range(4):
        gamma9half_over_sqrtpi *= Q(2*k+1, 2)
    integral0_coefficient = 2*factorial(1)
    integral3half_coefficient = 2*factorial(4)
    det_H_minus5_coefficient = Q(integral3half_coefficient, 32*factorial(4))
    U2_over_pi_coefficient = 64*det_H_minus5_coefficient*integral0_coefficient
    V_coefficient = 8*Q(integral3half_coefficient, 16)/gamma9half_over_sqrtpi
    assert gamma9half_over_sqrtpi == Q(105, 16)
    assert det_H_minus5_coefficient == Q(1, 16)
    assert U2_over_pi_coefficient == 8 and V_coefficient == Q(128, 35)
    # Budget square and the eigenvalue geometry; pi^(-3) is factored out.
    budget_square_coefficient = Q(8**3, 32**2)/Q(128, 35)**4
    assert budget_square_coefficient == 2*Q(1225, 32768)**2
    assert Q(16384**3, 4*12288**4) == Q(1, 144**2)
    assert Q(3, 2)*Q(-6)-2*Q(-9, 2) == 0
    point_upper = Q(1225, 32768)*Q(1, 3)*Q(5, 6)*Q(17, 15)*Q(8, 7)/144
    assert point_upper == Q(2975, 31850496) == Q(saved["point_budget_upper_rational"])
    assert point_upper < Q(1, 10000)
    epsilon, delta = Q(1, 1000), Q(1, 4000)
    eta = 9*epsilon
    assert 4**2*5 < 9**2 and eta < 1
    # The two profile values at x0 are k/12, k=1,5.
    for k in (1, 5):
        assert k >= 1 and k+7 <= 12
    factor = ((1+eta)/(1-eta))**9/(1-epsilon)**4
    assert factor == Q(saved["patch_factor"]) < Q(6, 5)
    patch_upper = factor*Q(1, 10000)
    assert patch_upper < Q(3, 25000) < Q(15198, 10**8)
    assert 4*delta**2 == Q(1, 2000)**2 < epsilon**2
    assert 5*delta**2 < epsilon**2
    # x0 is a half-sum of two allowed, non-collinear directions.
    v, w = (1, 1, 1, 1, 1), (1, 1, 1, -1, -1)
    x0 = tuple(Q(a+b, 2) for a, b in zip(v, w))
    assert x0 == (1, 1, 1, 0, 0) and x0 != v and x0 != w
    return {
        "direction_counts": direction_counts,
        "determinant_polynomial_low_to_high": list(map(str, determinant)),
        "bD_polynomial_low_to_high": list(map(str, bD)),
        "t_squared_lower_strict": "1/5", "t_squared_upper_strict": "11/7",
        "rho_symbolic_lower_difference": list(map(str, lower_difference)),
        "rho_symbolic_upper_difference": list(map(str, upper_difference)),
        "Mellin_U_squared_pi_coefficient": str(U2_over_pi_coefficient),
        "Mellin_V_coefficient": str(V_coefficient),
        "point_uniform_upper": str(point_upper),
        "relative_lambda_positive_rho_bound": "1/7",
        "relative_lambda_negative_rho_bound": "1",
        "epsilon": str(epsilon), "eta": str(eta), "patch_factor": str(factor),
        "strict_patch_cap": "3/25000", "threshold_compared": "15198/100000000",
        "explicit_face_patch_measure": str(delta**4),
        "explicit_cube_patch_measure": str(delta**5),
        "face_patch_max_distance_squared": str(4*delta**2),
        "cube_patch_max_distance_squared": str(5*delta**2),
        "x0_half_sum": [v, w],
        "no_rho_samples_or_old_fixed_rho_values_evaluated": True,
        "no_5_by_5_symbolic_matrix_recomputed": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--run", action="store_true", help="execute this one finite new audit; default hashes only")
    args = parser.parse_args()
    root = args.project_root.resolve()
    if (root/"cube-dim-5-bounds").is_dir():
        root = root/"cube-dim-5-bounds"
    regpath = HERE/"registration.json"
    reg = json.loads(regpath.read_text())
    assert sha(Path(__file__)) == reg["script_sha256"]
    sources = {}
    for source in reg["sources"]:
        path = root/source["project_relative"]
        assert sha(path) == source["sha256"], str(path)
        sources[source["id"]] = path
    if not args.run:
        print(json.dumps({"status": "HASH_ONLY_PASS", "sources": len(sources), "run": False}))
        return
    started = time.monotonic()
    print(json.dumps({"status": "AUDIT_STARTED", "actual_pid": os.getpid()}), flush=True)
    record = evaluate(sources)
    result = {"status": "PASS_CONTINUOUS_RHO_AND_POSITIVE_MEASURE_PATCH_BUDGET_OBSTRUCTION",
              "actual_pid": os.getpid(), "python": sys.executable,
              "arithmetic": "Python standard-library Fraction, integer polynomials and factorials; no SymPy/Arb",
              "registration_sha256": sha(regpath), "script_sha256": sha(Path(__file__)),
              "source_locks": reg["sources"], "record": record,
              "seconds": time.monotonic()-started,
              "true_G_upper_bound_claim": False, "new_C5_bound": False}
    (HERE/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "actual_pid": result["actual_pid"],
                      "seconds": result["seconds"], "results_sha256": sha(HERE/"results.json")}))


if __name__ == "__main__":
    main()
