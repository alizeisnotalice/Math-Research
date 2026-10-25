"""Three conditional Young minorants, exact mixed moments, rational Bernstein proof."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
from math import comb, factorial, floor, ceil
import hashlib
import importlib.util
import json
import sys
import time
import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'astra_upper001' / 'verify.py'
spec = importlib.util.spec_from_file_location('upper001_locked', SOURCE)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
BETA = F(3, 4)
ORDER = 11
PLUS = {i for i, v in enumerate(old.V) if v[0] * v[1] == 1}


def ar(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def smul(a, b):
    out = [F(0)] * (ORDER + 1)
    for i, ai in enumerate(a):
        if ai:
            for j in range(min(len(b), ORDER + 1 - i)):
                out[i + j] += ai * b[j]
    return out


@lru_cache(None)
def linear_power(exponent, sign):
    return tuple(old.choose(exponent, n) * sign ** n for n in range(ORDER + 1))


@lru_cache(None)
def denominator_power(exponent):
    # h(z)=((1-z)^beta+(1+z)^beta)/2, h(0)=1.
    delta = [F(0)] + [old.choose(BETA, n) if n % 2 == 0 else F(0) for n in range(1, ORDER + 1)]
    out = [F(0)] * (ORDER + 1)
    power = [F(1)] + [F(0)] * ORDER
    for j in range(ORDER // 2 + 1):
        factor = old.choose(F(exponent), j)
        out = [a + factor * b for a, b in zip(out, power)]
        power = smul(power, delta)
    return tuple(out)


@lru_cache(None)
def mixed_weight_moments(n1p, n2p, n1m, n2m):
    """E_R[prod P_i^k_i * r^n], r=sum_plus P - sum_minus P, 0<=n<=11."""
    kp, km = n1p + 2 * n2p, n1m + 2 * n2m
    degree = kp + km
    accum = [F(0)] * (ORDER + 1)
    for lp in range(n2p + 1):
        for lm in range(n2m + 1):
            jp = n1p + n2p + lp
            jm = n1m + n2m + lm
            j = jp + jm
            factor = (F(comb(n2p, lp) * comb(n2m, lm) * factorial(j), 16 ** j)
                      * BETA ** (n1p + n1m + n2p + n2m + lp + lm)
                      * (1 - BETA) ** (n2p + n2m - lp - lm)
                      / old.rising(BETA, degree))
            power = smul(linear_power(BETA * jp - kp, -1), linear_power(BETA * jm - km, 1))
            power = smul(power, denominator_power(-j - 1))
            accum = [a + factor * b for a, b in zip(accum, power)]
    result = tuple(accum[n] * factorial(n) / old.rising(BETA + degree, n) for n in range(ORDER + 1))
    expected_zero = old.laplace_moment([1] * (n1p + n1m) + [2] * (n2p + n2m))
    assert result[0] == expected_zero
    return result


def determinant_mixed_moments():
    subsets = {}
    first_coefficients = [0] * 6
    nonzero = 0
    for basis in combinations(range(16), 5):
        weight = old.bareiss_det([old.V[i] for i in basis]) ** 2
        if not weight:
            continue
        nonzero += 1
        p = len(set(basis) & PLUS)
        first_coefficients[p] += weight
        for size in range(6):
            for subset in combinations(basis, size):
                if subset not in subsets:
                    subsets[subset] = [0] * 6
                subsets[subset][p] += weight
    assert nonzero == 3008
    assert sum(first_coefficients) == 16 ** 5
    # b[jp,jm,l]=sum_{J with these signs} sum_{B,C contain J,p(B)+p(C)=l}wB*wC.
    b = {(jp, jm): [0] * 11 for jp in range(6) for jm in range(6 - jp)}
    for subset, coef in subsets.items():
        jp = len(set(subset) & PLUS)
        jm = len(subset) - jp
        for p, wp in enumerate(coef):
            if wp:
                for q, wq in enumerate(coef):
                    if wq:
                        b[jp, jm][p + q] += wp * wq
    patterns = {}
    for kp in range(6):
        for km in range(6 - kp):
            weights = [0] * 11
            for jp in range(kp, 6):
                for jm in range(km, 6 - jp):
                    sign = (-1) ** (jp + jm - kp - km) * comb(jp, kp) * comb(jm, km)
                    weights = [w + sign * z for w, z in zip(weights, b[jp, jm])]
            assert all(w >= 0 for w in weights)
            for total_plus_degree, weight in enumerate(weights):
                if weight:
                    pattern = (total_plus_degree - 2 * kp, kp,
                               10 - total_plus_degree - 2 * km, km)
                    assert min(pattern) >= 0
                    patterns[pattern] = patterns.get(pattern, 0) + weight
    assert sum(patterns.values()) == 16 ** 10
    first = [F(0)] * (ORDER + 1)
    second = [F(0)] * (ORDER + 1)
    for p, weight in enumerate(first_coefficients):
        if weight:
            first = [x + weight * y for x, y in zip(first, mixed_weight_moments(p, 0, 5 - p, 0))]
    for pattern, weight in patterns.items():
        second = [x + weight * y for x, y in zip(second, mixed_weight_moments(*pattern))]
    assert first[0] == F(648, 1463)
    assert second[0] == F(126404864, 474619145)
    assert all(first[n] == second[n] == 0 for n in range(1, ORDER + 1, 2))
    overlap_totals = [sum(w for pattern, w in patterns.items() if pattern[1] + pattern[3] == k) for k in range(6)]
    expected = [126508597248, 408315494400, 403408158720, 142857994240, 17825792000, 595591168]
    assert overlap_totals == expected
    return first, second, {'nonzero_bases': nonzero,
                           'basis_coefficients_by_plus_count': first_coefficients,
                           'pair_pattern_weights': {','.join(map(str, k)): v for k, v in sorted(patterns.items())},
                           'overlap_totals_match_001': True}


def cofactor_mixed_moments():
    def cb(n, k):
        return comb(n, k) if 0 <= k <= n else 0
    # A(a,b)=sum_plus a vv'+sum_minus b vv'; diag 8(a+b), A12=8(a-b).
    coefficients = {
        'B11': [8 ** 4 * cb(4, p) for p in range(5)],
        'B33': [4 * 8 ** 4 * cb(2, p - 1) for p in range(5)],
        'B12': [-8 ** 4 * (cb(3, p - 1) - cb(3, p)) for p in range(5)],
    }
    out = {}
    for name, coef in coefficients.items():
        row = [F(0)] * (ORDER + 1)
        for p, weight in enumerate(coef):
            row = [x + weight * y for x, y in zip(row, mixed_weight_moments(p, 0, 4 - p, 0))]
        out[name] = row
    assert out['B11'][0] == out['B33'][0] == F(216, 385)
    assert out['B12'][0] == 0
    assert all(out['B11'][n] == out['B33'][n] == 0 for n in range(1, ORDER + 1, 2))
    assert all(out['B12'][n] == 0 for n in range(0, ORDER + 1, 2))
    assert all(out[name][n] >= 0 for name in ['B11', 'B33'] for n in range(0, ORDER + 1, 2))
    return out, coefficients


def bernstein_coefficients(power, left, right):
    degree = len(power) - 1
    local = [sum((power[i] * comb(i, j) * left ** (i - j) for i in range(j, degree + 1)), F(0))
             * (right - left) ** j for j in range(degree + 1)]
    return [sum((local[j] * F(comb(k, j), comb(degree, j)) for j in range(k + 1)), F(0)) for k in range(degree + 1)]


def certify_lower_bound(power, target, max_depth=16):
    residual = list(power)
    residual[0] -= target
    stack = [(F(0), F(1), 0)]
    leaves = []
    while stack:
        left, right, depth = stack.pop()
        coeff = bernstein_coefficients(residual, left, right)
        if min(coeff) >= 0:
            leaves.append({'left': str(left), 'right': str(right), 'depth': depth,
                           'minimum_residual_Bernstein_coefficient': str(min(coeff))})
        elif depth >= max_depth:
            return False, leaves, {'left': str(left), 'right': str(right), 'depth': depth}
        else:
            mid = (left + right) / 2
            stack.extend([(mid, right, depth + 1), (left, mid, depth + 1)])
    ordered = sorted(leaves, key=lambda row: F(row['left']))
    assert ordered[0]['left'] == '0' and ordered[-1]['right'] == '1'
    assert all(F(a['right']) == F(b['left']) for a, b in zip(ordered, ordered[1:]))
    return True, ordered, None


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def ppow(a, n):
    out = [F(1)]
    for _ in range(n):
        out = pmul(out, a)
    return out


def bernstein_family_exclusion(U, V, c_upper):
    v7, u11 = ppow(V, 7), ppow(U, 11)
    degree = max(len(v7), len(u11)) - 1
    residual = [F(0)] * (degree + 1)
    for i, value in enumerate(v7):
        residual[i] += F(1, 27000) ** 4 * value
    for i, value in enumerate(u11):
        residual[i] -= c_upper ** 4 * value
    initial = bernstein_coefficients(residual, F(0), F(1))
    stack = [(F(0), F(1), 0, initial)]
    leaves = []
    while stack:
        left, right, depth, coeff = stack.pop()
        if min(coeff) > 0:
            leaves.append({'left': str(left), 'right': str(right), 'depth': depth,
                           'minimum_Bernstein_coefficient': str(min(coeff))})
            continue
        if depth == 10:
            return {'certified': False, 'degree': degree, 'failed_interval': [str(left), str(right)]}
        work = coeff[:]
        left_coeff, right_coeff = [work[0]], [work[-1]]
        for _ in range(degree):
            work = [(a + b) / 2 for a, b in zip(work, work[1:])]
            left_coeff.append(work[0])
            right_coeff.append(work[-1])
        right_coeff.reverse()
        mid = (left + right) / 2
        stack.extend([(mid, right, depth + 1, right_coeff), (left, mid, depth + 1, left_coeff)])
    leaves.sort(key=lambda row: F(row['left']))
    assert leaves[0]['left'] == '0' and leaves[-1]['right'] == '1'
    assert all(a['right'] == b['left'] for a, b in zip(leaves, leaves[1:]))
    return {'certified': True, 'degree': degree, 'target_ceiling': '1/27000',
            'parameter_domain': 'all a>=0 and all |eta|<=1', 'Bernstein_partition': leaves,
            'residual_power_coefficients': [str(x) for x in residual]}


def main():
    began = time.monotonic()
    ctx.prec = 256
    d1, d2, geometry = determinant_mixed_moments()
    cofactor, cofactor_coeff = cofactor_mixed_moments()
    s = F(5, 11)
    tau = d2[0] / d1[0]
    A = ar(2 - s) * ar(tau) ** ar(s - 1)
    B = ar(1 - s) * ar(tau) ** ar(s - 2)
    c = 3 / (2 ** (arb(15) / 2) * arb.pi() ** (arb(5) / 2))
    A_lower = F(floor(float(A) * 10 ** 12) - 1, 10 ** 12)
    B_upper = F(ceil(float(B) * 10 ** 12) + 1, 10 ** 12)
    c_lower = F(floor(float(c) * 10 ** 15) - 1, 10 ** 15)
    A_upper = F(ceil(float(A) * 10 ** 12) + 1, 10 ** 12)
    B_lower = F(floor(float(B) * 10 ** 12) - 1, 10 ** 12)
    c_upper = F(ceil(float(c) * 10 ** 15) + 1, 10 ** 15)
    assert ar(A_lower) < A and ar(B_upper) > B and ar(c_lower) < c
    assert ar(A_upper) > A and ar(B_lower) < B and ar(c_upper) > c
    assert A_lower > B_upper > 0
    ell = [A_lower * x - B_upper * y for x, y in zip(d1, d2)]
    assert all(ell[n] >= 0 for n in range(0, 8, 2))
    a = F(2, 3)
    candidates = []
    for eta in [F(1, 8), F(1, 4), F(1, 2)]:
        power = [F(0)] * 7
        for k in range(0, 8, 2):
            power[k // 2] += F(11, 4) * a ** 7 * comb(7, k) * eta ** k * ell[k]
        for k in range(0, 12, 2):
            common = -F(7, 4) * a ** 11 * comb(11, k) * eta ** k
            power[k // 2] += common * (cofactor['B11'][k] + 3 * cofactor['B33'][k])
            power[k // 2 + 1] += common * cofactor['B11'][k]
        for k in range(1, 12, 2):
            power[(k + 1) // 2] -= F(7, 2) * a ** 11 * comb(11, k) * eta ** k * cofactor['B12'][k]
        approximate_roots = np.polynomial.polynomial.polyroots([float(k * power[k]) for k in range(1, 7)])
        points = [0.0, 1.0] + [float(z.real) for z in approximate_roots if abs(z.imag) < 1e-8 and 0 <= z.real <= 1]
        evaluations = [float(np.polynomial.polynomial.polyval(x, [float(v) for v in power])) for x in points]
        proxy_min = min(evaluations)
        target = F(floor(proxy_min * float(c_lower) * 10 ** 8) - 1, 10 ** 8)
        row = {'eta': str(eta), 'a': str(a), 'global_space_polynomial_power_coefficients': [str(x) for x in power],
               'float_min_only_for_target_selection': proxy_min,
               'float_argmin_only_diagnostic': points[evaluations.index(proxy_min)],
               'float_g_proxy_not_certificate': proxy_min * float(c_lower)}
        corner_U = sum((F(comb(7, k)) * eta ** k * ell[k] for k in range(0, 8, 2)), F(0))
        corner_V = sum((F(comb(11, k)) * eta ** k * (2 * cofactor['B11'][k] + 3 * cofactor['B33'][k]) for k in range(0, 12, 2)), F(0))
        corner_V += sum((2 * F(comb(11, k)) * eta ** k * cofactor['B12'][k] for k in range(1, 12, 2)), F(0))
        assert corner_U > 0 and corner_V > 0
        ceiling = c * ar(corner_U) ** (arb(11) / 4) / ar(corner_V) ** (arb(7) / 4)
        row.update({'exact_corner_U': str(corner_U), 'exact_corner_V': str(corner_V),
                    'all_a_corner_ceiling_Arb': str(ceiling),
                    'all_a_ceiling_below_1_over_27000': bool(ceiling < arb(1) / 27000),
                    'nonconstant_space_coefficients_all_negative': all(x < 0 for x in power[1:])})
        if target > 0:
            ok, leaves, failure = certify_lower_bound(power, target / c_lower)
            attempts = 1
            if not ok:
                target -= F(1, 10 ** 8)
                ok, leaves, failure = certify_lower_bound(power, target / c_lower)
                attempts += 1
            row.update({'certified': ok, 'target_g': str(target), 'attempts': attempts,
                        'strictly_improves_1_over_27000': ok and target > F(1, 27000),
                        'Bernstein_partition': leaves, 'failed_interval': failure})
        else:
            # A finite exact witness that this polynomial cannot certify a positive bound.
            witness = next((F(k, 64) for k in range(65) if sum(power[j] * F(k, 64) ** j for j in range(7)) < 0), None)
            row.update({'certified': False, 'positive_target_rejected': True,
                        'negative_rational_witness': str(witness) if witness is not None else None})
        candidates.append(row)
    U_upper = [F(0)] * 7
    V_exact = [F(0)] * 12
    for k in range(0, 8, 2):
        U_upper[k] = comb(7, k) * (A_upper * d1[k] - B_lower * d2[k])
    for k in range(0, 12, 2):
        V_exact[k] = comb(11, k) * (2 * cofactor['B11'][k] + 3 * cofactor['B33'][k])
    for k in range(1, 12, 2):
        V_exact[k] = 2 * comb(11, k) * cofactor['B12'][k]
        assert V_exact[k] <= 0
    family_exclusion = bernstein_family_exclusion(U_upper, V_exact, c_upper)
    out = {'status': 'three_registered_conditional_Young_polynomial_candidates', 'new_C5_upper_bound': False,
           'constants': {'A_arb': str(A), 'B_arb': str(B), 'A_lower': str(A_lower), 'B_upper': str(B_upper), 'c_lower': str(c_lower),
                         'A_upper': str(A_upper), 'B_lower': str(B_lower), 'c_upper': str(c_upper)},
           'geometry': geometry, 'cofactor_polynomial_coefficients': cofactor_coeff,
           'exact_mixed_moments': {'E_d_r_power': [str(x) for x in d1], 'E_d_squared_r_power': [str(x) for x in d2],
                                   **{f'E_{name}_r_power': [str(x) for x in row] for name, row in cofactor.items()}},
           'candidates': candidates,
           'whole_linear_single_correlation_family_exclusion': family_exclusion,
           'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE, ROOT / 'registration.md', Path(__file__)]},
           'runtime': {'python': sys.executable, 'arb_bits': ctx.prec, 'seconds': time.monotonic() - began}}
    (ROOT / 'verification.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'candidates': candidates, 'family_exclusion_summary': {k: v for k, v in family_exclusion.items() if k not in ['residual_power_coefficients','Bernstein_partition']},
                      'family_exclusion_Bernstein_leaves': len(family_exclusion.get('Bernstein_partition', [])), 'runtime': out['runtime']}, indent=2))


if __name__ == '__main__':
    main()
