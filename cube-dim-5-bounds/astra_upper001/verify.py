"""Finite exact checks and Arb enclosure; no sampling and no C5 claim."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product, combinations
from math import factorial, comb
import hashlib
import json
import platform
import sys
import time

import sympy as sp
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5')
BETA = F(3, 4)
V = [(1,) + s for s in product((-1, 1), repeat=4)]


def rising(a, k):
    out = F(1)
    for j in range(k):
        out *= a + j
    return out


def choose(a, k):
    out = F(1)
    for j in range(k):
        out *= (a - j) / (j + 1)
    return out


def transformed_moment(columns, counts):
    """Exact coefficient of 16/sum_i(1-sum_j a_ij*z_j)^beta."""
    shape = tuple(counts)
    indices = sorted(product(*(range(k + 1) for k in shape)), key=sum)
    zero = (0,) * len(shape)
    den = {zero: F(1)}
    inv = {zero: F(1)}
    for powers in indices[1:]:
        degree = sum(powers)
        multinomial = F(factorial(degree))
        for power in powers:
            multinomial /= factorial(power)
        avg = F(0)
        for i in range(16):
            term = F(1)
            for col, power in zip(columns, powers):
                term *= F(col[i]) ** power
            avg += term / 16
        den[powers] = choose(BETA, degree) * (-1) ** degree * multinomial * avg
        subtotal = F(0)
        for part in product(*(range(k + 1) for k in powers)):
            if part != zero:
                rest = tuple(k - j for k, j in zip(powers, part))
                subtotal += den[part] * inv[rest]
        inv[powers] = -subtotal
    out = inv[shape] / rising(BETA, sum(shape))
    for k in shape:
        out *= factorial(k)
    return out


def laplace_coefficients(k):
    """(-d/dt)^k exp(-t^beta)=exp(-t^beta) sum_j a[k,j]t^(j beta-k)."""
    coef = {0: F(1)}
    for order in range(k):
        nxt = {}
        for j, value in coef.items():
            nxt[j] = nxt.get(j, F(0)) + (order - j * BETA) * value
            nxt[j + 1] = nxt.get(j + 1, F(0)) + BETA * value
        coef = {j: val for j, val in nxt.items() if val}
    return coef


def laplace_moment(counts):
    coef = {0: F(1)}
    for k in counts:
        nxt = {}
        for j1, a in coef.items():
            for j2, b in laplace_coefficients(k).items():
                nxt[j1 + j2] = nxt.get(j1 + j2, F(0)) + a * b
        coef = nxt
    return sum((a * factorial(j) / 16 ** j for j, a in coef.items()), F(0)) / rising(BETA, sum(counts))


def bareiss_det(rows):
    a = [list(row) for row in rows]
    sign, old = 1, 1
    for k in range(len(a) - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, len(a)) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        current = a[k][k]
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                value = a[i][j] * current - a[i][k] * a[k][j]
                assert value % old == 0
                a[i][j] = value // old
            a[i][k] = 0
        old = current
    return sign * a[-1][-1]


def determinant_squared_moment():
    # det M = sum_B det(V_B)^2 prod_{i in B}P_i (Cauchy--Binet).
    # b_j=sum_J (sum_{B contains J}w_B)^2 = sum_k choose(k,j)a_k.
    supersets = [{} for _ in range(6)]
    nonzero, total_weight = 0, 0
    for basis in combinations(range(16), 5):
        determinant = bareiss_det([V[i] for i in basis])
        weight = determinant ** 2
        if not weight:
            continue
        nonzero += 1
        total_weight += weight
        for j in range(6):
            for subset in combinations(basis, j):
                supersets[j][subset] = supersets[j].get(subset, 0) + weight
    assert total_weight == 16 ** 5
    b = [sum(w ** 2 for w in s.values()) for s in supersets]
    a = [sum((-1) ** (j - k) * comb(j, k) * b[j] for j in range(k, 6)) for k in range(6)]
    assert all(x >= 0 for x in a)
    assert sum(a) == total_weight ** 2
    overlap_moments = [laplace_moment([2] * k + [1] * (10 - 2 * k)) for k in range(6)]
    result = sum((weight * moment for weight, moment in zip(a, overlap_moments)), F(0))
    return {'enumerated_subsets': comb(16, 5), 'nonzero_bases': nonzero,
            'cauchy_binet_coefficient_sum': total_weight,
            'ordered_pair_overlap_weight_sums': a,
            'overlap_moments': [str(x) for x in overlap_moments],
            'E_detM_squared': str(result)}, result


def main():
    began = time.monotonic()
    ctx.prec = 256
    frame = sp.Matrix(V).T * sp.Matrix(V)
    assert frame == 16 * sp.eye(5)
    basis = [
        (1, 1, 1, 1, 1),
        (1, -1, -1, -1, 1),
        (1, -1, -1, 1, -1),
        (1, -1, 1, -1, -1),
        (1, 1, -1, -1, -1),
    ]
    assert abs(sp.Matrix(basis).det()) == 48
    covers = {v: 0 for v in V}
    for signs in V:
        reflected = [tuple(a * b for a, b in zip(signs, v)) for v in basis]
        assert abs(sp.Matrix(reflected).det()) == 48
        for v in reflected:
            covers[v] += 1
    assert set(covers.values()) == {5}

    pi = arb.pi()
    one_d_peak = (arb(2) / 3).gamma() * 2 / (3 * pi)
    peak_ceiling = one_d_peak ** 5 / 48 * (arb(5) / 16) ** (arb(10) / 3)
    kappa = 3 * (arb(7) / 4).gamma() / (2 ** (arb(7) / 2) * pi ** (arb(5) / 2))
    ed_ceiling = (4 * pi) ** (arb(5) / 2) * peak_ceiling
    trace_ceiling = kappa * 5 ** (-arb(7) / 4) * ed_ceiling ** (arb(3) / 10)
    assert trace_ceiling < arb(84) / 1000000 < arb(1) / 10000
    normalization = 1 / (16 * (arb(7) / 4).gamma())
    normalized_prefactor = 3 / (2 ** (arb(15) / 2) * pi ** (arb(5) / 2))
    assert (kappa * normalization - normalized_prefactor).contains(0)

    det_receipt, expected_det_squared = determinant_squared_moment()
    spectral_ceiling_fourth = F(5) ** 7 * F(15, 43) ** 3 * F(50, 43) ** 40
    assert spectral_ceiling_fourth < 35 ** 4
    det_minorant_lower = normalized_prefactor * arb(expected_det_squared.numerator) / expected_det_squared.denominator / 35

    chars = {(i, j): [v[i] * v[j] for v in V] for i in range(5) for j in range(i + 1, 5)}
    covariance = [
        [transformed_moment([a, b], [1, 1]) for b in chars.values()]
        for a in chars.values()
    ]
    assert covariance == [[F(int(i == j), 7) for j in range(10)] for i in range(10)]
    single = [transformed_moment([chars[(0, 1)]], [k]) for k in range(7)]
    joint = {
        'E_r12_r13_r23': transformed_moment([chars[(0, 1)], chars[(0, 2)], chars[(1, 2)]], [1, 1, 1]),
        'E_r12_squared_r13_squared': transformed_moment([chars[(0, 1)], chars[(0, 2)]], [2, 2]),
        'E_r12_squared_r34_squared': transformed_moment([chars[(0, 1)], chars[(2, 3)]], [2, 2]),
    }
    checks = {}
    for counts in ([1], [2], [3], [4], [1, 1], [2, 1], [3, 1], [2, 2], [1, 1, 1, 1]):
        indicators = [[int(i == j) for i in range(16)] for j in range(len(counts))]
        from_transform = transformed_moment(indicators, counts)
        from_laplace = laplace_moment(counts)
        assert from_transform == from_laplace
        checks[','.join(map(str, counts))] = str(from_transform)

    source_names = [
        'AGENTS.md',
        'cube-dim-5-bounds/alternative_alpha/transfer.md',
        'cube-dim-5-bounds/upper/round009_alpha_shorttime.md',
        'cube-dim-5-bounds/round005/alpha15_middle_budget.md',
        'cube-dim-5-bounds/round004/gaussian_mixture_lower.md',
        'cube-dim-5-bounds/round004/isotropic_exclusion.md',
    ]
    sources = {str(MAIN / name): hashlib.sha256((MAIN / name).read_bytes()).hexdigest() for name in source_names}
    sources[str(ROOT / 'verify.py')] = hashlib.sha256((ROOT / 'verify.py').read_bytes()).hexdigest()
    sources[str(ROOT / 'registration.md')] = hashlib.sha256((ROOT / 'registration.md').read_bytes()).hexdigest()
    sources[str(ROOT / 'report.md')] = hashlib.sha256((ROOT / 'report.md').read_bytes()).hexdigest()
    out = {
        'status': 'rigorous_trace_mechanism_exclusion_and_exact_moment_interface_checks',
        'C5_upper_bound_improved': False,
        'g_target_achieved': False,
        'exact_geometry': {
            'frame': [[int(x) for x in row] for row in frame.tolist()],
            'basis_rows': basis, 'basis_determinant': int(sp.Matrix(basis).det()),
            'reflected_bases': 16, 'each_direction_occurrences': 5,
        },
        'arb_256bit': {
            'one_dimensional_peak': str(one_d_peak),
            'p1_zero_upper_bound': str(peak_ceiling),
            'kappa': str(kappa),
            'trace_jensen_mechanism_ceiling': str(trace_ceiling),
            'certified_ceiling_lt': '84/1000000',
            'g_target': '1/10000',
            'ET_power_minus_beta': str(normalization),
            'normalized_shape_prefactor': str(normalized_prefactor),
            'det_squared_minorant_global_g': str(det_minorant_lower),
        },
        'finite_determinant_minorant': {
            **det_receipt,
            'minorant': 'H(M,x)>=det(M)^2/35 for all x in cube and M positive definite with trace 5',
            'exact_spectral_fourth_power_ceiling': str(spectral_ceiling_fourth),
            'spectral_ceiling_strictly_below_35_fourth': True,
            'status': 'certified_global_but_weaker_than_existing_tilted_bound',
        },
        'exact_moments': {
            'single_correlation_orders_0_through_6': [str(x) for x in single],
            'ten_correlation_covariance': [[str(x) for x in row] for row in covariance],
            'E_frobenius_M_minus_identity_squared': str(F(20, 7)),
            'selected_joint_moments': {key: str(val) for key, val in joint.items()},
            'independent_laplace_weight_moment_crosschecks': checks,
        },
        'runtime': {'python': sys.executable, 'version': platform.python_version(), 'sympy': sp.__version__, 'arb_precision_bits': ctx.prec},
        'source_sha256': sources,
        'elapsed_seconds': time.monotonic() - began,
    }
    (ROOT / 'verification.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
