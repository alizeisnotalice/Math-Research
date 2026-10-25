#!/usr/bin/env python3
"""Exact finite algebra supporting the analytic theorem; zero kernel passes."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from time import perf_counter
import hashlib
import json
import math
import os
import signal
import sys
import flint
from flint import fmpq, fmpq_mat

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
SOURCES = [
    PROJECT / 'upper/round006/g/sharp_integrand_ode.md',
    PROJECT / 'upper/round008/chamber_convexity/report.md',
    Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round009/short_potential/curvature_results.json'),
    Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper/round009/short_potential_barrier/theorem.md'),
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n')


def timeout(signum, frame):
    raise TimeoutError('registered 1200-second structural cap')


def main():
    start = perf_counter()
    c = [Q(1)] + [Q(1, 2)] * 5
    z = [Q(n, 200) for n in [0, 51, 22, 2, -24, -51]]
    inputs = [HERE / 'registration.md', HERE / 'theorem.md', Path(__file__)] + SOURCES
    manifest = {
        'status': 'REGISTERED_BEFORE_EXACT_STRUCTURAL_VERIFICATION',
        'T_of_interest': '4/125', 'theorem_T_scope': 'every real T>0',
        'c': list(map(str, c)), 'z': list(map(str, z)),
        'max_seconds': 1200, 'max_new_H_passes': 0,
        'max_new_G_passes': 0, 'max_new_p_passes': 0,
        'inputs': [{'path': str(p), 'sha256': sha(p)} for p in inputs],
    }
    write('input_manifest.json', manifest)
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(1200)
    try:
        V = [(1,) + s for s in product([-1, 1], repeat=5)]
        A = fmpq_mat([[v[j] for v in V] for j in range(6)])
        gram = A * A.transpose()
        assert all(gram[i, j] == (32 if i == j else 0)
                   for i in range(6) for j in range(6))
        v = (1,) * 6
        w = (1,) + (-1,) * 5
        basis = [v, w]
        for j in range(1, 5):
            vv = list(v)
            vv[j] = -1
            basis.append(tuple(vv))
        assert len(set(basis)) == 6 and all(b in V for b in basis)
        B = fmpq_mat([[b[j] for b in basis] for j in range(6)])
        det = B.det()
        assert det != 0
        inv = B.inv()
        cc = fmpq_mat([[fmpq(x.numerator, x.denominator)] for x in c])
        solved = inv * cc
        expected = [Q(3, 4), Q(1, 4), Q(0), Q(0), Q(0), Q(0)]
        assert all(solved[j, 0] == fmpq(q.numerator, q.denominator)
                   for j, q in enumerate(expected))
        rest = [col for col in V if col not in basis]
        assert len(rest) == 26
        N = fmpq_mat([[col[j] for col in rest] for j in range(6)])
        M = inv * N
        row_bounds = [sum((abs(M[i, j]) for j in range(26)), fmpq(0))
                      for i in range(6)]
        C = max(row_bounds)
        assert C > 0
        t0 = min(fmpq(1), 1 / C)
        lower_times_pi6 = fmpq(1) / (2 ** 26 * 25 * abs(det) * (1 + C ** 2) ** 4)
        assert lower_times_pi6 > 0 and t0 > 0

        endpoint = [a + b for a, b in zip(c, z)]
        assert endpoint[0] == 1 and endpoint[-1] > 0
        assert all(endpoint[j] > endpoint[j + 1] for j in range(5))
        assert sum(z[1:]) == 0
        assert all(z[j] > z[j + 1] for j in range(1, 5))
        partials = [sum(z[1:k + 1]) for k in range(1, 5)]
        assert partials == [Q(n, 200) for n in [51, 73, 75, 51]]
        # Linear coefficient check of the exact summation-by-parts identity.
        rebuilt = [partials[0]] + [partials[j] - partials[j - 1]
                                    for j in range(1, 4)] + [-partials[-1]]
        assert rebuilt == z[1:]
        walls = []
        for i in range(1, 5):
            m = endpoint.copy()
            m[i] = m[i + 1] = (m[i] + m[i + 1]) / 2
            assert len(set(m)) == 5 and min(m) > 0
            walls.append({'one_based_pair': [i + 1, i + 2],
                          'endpoint_wall': list(map(str, m)),
                          'distinct_absolute_values': len(set(map(abs, m)))})
        patterns = list(product([-1, 1], repeat=3))
        representatives = {s if s[0] == 1 else tuple(-a for a in s)
                           for s in patterns}
        assert len(patterns) == 8 and len(representatives) == 4

        # r=(u/t)^2 >= 0: |6r-2| <= 2(1+r)^2.
        # Upper difference = 2(r-1/2)^2+7/2; lower difference = 2r^2+10r.
        upper_coefficients = [Q(4), Q(-2), Q(2)]
        completed_square = [2 * Q(1, 2) ** 2 + Q(7, 2), -4 * Q(1, 2), Q(2)]
        assert upper_coefficients == completed_square
        lower_coefficients = [Q(0), Q(10), Q(2)]
        assert all(a >= 0 for a in lower_coefficients)

        out = {
            'status': 'PASS_EXACT_STRUCTURAL_CHECK_ONLY',
            'pid': os.getpid(), 'python': sys.version,
            'python_executable': sys.executable, 'python_flint': flint.__version__,
            'arithmetic': 'Fraction and FLINT exact rational matrices; no floating point kernel',
            'V_cardinality': len(V), 'AA_transpose_is_32_identity': True,
            'basis_columns': basis, 'basis_determinant': str(det),
            'B_inverse_c': [str(solved[j, 0]) for j in range(6)],
            'remaining_columns': len(rest), 'row_l1_bounds_B_inverse_A_N': list(map(str, row_bounds)),
            'C': str(C), 't0': str(t0),
            'positive_p_t_c_lower_coefficient_times_pi6': str(lower_times_pi6),
            'density_lower_power': 't^-2',
            'c': list(map(str, c)), 'z': list(map(str, z)),
            'x_1': list(map(str, endpoint)),
            'all_0_lt_s_le_1_path_in_strict_chamber': True,
            'z_sum_last_five': str(sum(z[1:])), 'z_positive_partial_sums': list(map(str, partials)),
            'summation_by_parts_coefficients_verified': True,
            'single_equality_wall_examples': walls,
            'three_sign_patterns': patterns, 'opposite_pattern_representatives': sorted(representatives),
            'S3_max_distinct_absolute_values': 4,
            'second_Cauchy_relative_derivative_bound': 2,
            'bound_polynomial_upper_difference': list(map(str, upper_coefficients)),
            'bound_polynomial_lower_difference': list(map(str, lower_coefficients)),
            'new_H_passes': 0, 'new_G_passes': 0, 'new_p_passes': 0,
            'numeric_negative_Jensen_triple_located': False,
            'analytic_nonconvexity_proof': 'theorem.md',
            'global_m_T_certified': False, 'new_C6_bound': False,
            'source_sha256': sha(__file__),
            'manifest_sha256': sha(HERE / 'input_manifest.json'),
            'elapsed_seconds': perf_counter() - start,
        }
        write('structural_certificate.json', out)
        print(json.dumps({k: out[k] for k in ['status', 'basis_determinant', 'C', 't0',
            'positive_p_t_c_lower_coefficient_times_pi6', 'new_H_passes', 'new_G_passes',
            'new_p_passes', 'elapsed_seconds']}, indent=2))
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
