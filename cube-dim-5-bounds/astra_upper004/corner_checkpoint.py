"""Necessary corner checks only. No global g claim from these evaluations."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, product
from functools import lru_cache
from math import comb, factorial, floor
import hashlib
import importlib.util
import json
import sys
import time
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'astra_upper001' / 'verify.py'
spec = importlib.util.spec_from_file_location('upper001_locked', SOURCE)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
BETA, ORDER = F(3, 4), 11
WEIGHTS, COUNTS = (10, 2, -2), (1, 5, 10)
GROUP = [WEIGHTS.index((sum(v) ** 2 - 5) // 2) for v in old.V]
assert tuple(GROUP.count(i) for i in range(3)) == COUNTS


def ar(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def smul(a, b):
    out = [F(0)] * (ORDER + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(min(len(b), ORDER + 1 - i)):
                out[i + j] += x * b[j]
    return out


@lru_cache(None)
def linear_power(exponent, weight):
    return tuple(old.choose(exponent, n) * (-weight) ** n for n in range(ORDER + 1))


@lru_cache(None)
def denominator_power(exponent):
    delta = [F(0)] + [sum((F(count, 16) * old.choose(BETA, n) * (-weight) ** n for weight, count in zip(WEIGHTS, COUNTS)), F(0)) for n in range(1, ORDER + 1)]
    assert delta[1] == 0
    out = [F(0)] * (ORDER + 1)
    power = [F(1)] + [F(0)] * ORDER
    for j in range(ORDER // 2 + 1):
        out = [x + old.choose(F(exponent), j) * y for x, y in zip(out, power)]
        power = smul(power, delta)
    return tuple(out)


@lru_cache(None)
def mixed(n1, n2):
    degrees = tuple(a + 2 * b for a, b in zip(n1, n2))
    degree = sum(degrees)
    out = [F(0)] * (ORDER + 1)
    for ell in product(*(range(k + 1) for k in n2)):
        js = tuple(a + b + e for a, b, e in zip(n1, n2, ell))
        j = sum(js)
        factor = F(factorial(j), 16 ** j) / old.rising(BETA, degree)
        for k, e in zip(n2, ell):
            factor *= comb(k, e)
        factor *= BETA ** (sum(n1) + sum(n2) + sum(ell)) * (1 - BETA) ** (sum(n2) - sum(ell))
        power = list(denominator_power(-j - 1))
        for count_j, count_k, weight in zip(js, degrees, WEIGHTS):
            power = smul(power, linear_power(BETA * count_j - count_k, weight))
        out = [x + factor * y for x, y in zip(out, power)]
    moments = tuple(out[n] * factorial(n) / old.rising(BETA + degree, n) for n in range(ORDER + 1))
    assert moments[0] == old.laplace_moment([1] * sum(n1) + [2] * sum(n2))
    return moments


def counts(indices):
    return tuple(sum(GROUP[i] == group for i in indices) for group in range(3))


def determinant_moments():
    containing, first_coef = {}, {}
    for basis in combinations(range(16), 5):
        weight = old.bareiss_det([old.V[i] for i in basis]) ** 2
        if not weight:
            continue
        profile = counts(basis)
        first_coef[profile] = first_coef.get(profile, 0) + weight
        for size in range(6):
            for subset in combinations(basis, size):
                group = containing.setdefault(subset, {})
                group[profile] = group.get(profile, 0) + weight
    overlap_binomial = {}
    for subset, coef in containing.items():
        profile = counts(subset)
        dest = overlap_binomial.setdefault(profile, {})
        for a, wa in coef.items():
            for b, wb in coef.items():
                total = tuple(x + y for x, y in zip(a, b))
                dest[total] = dest.get(total, 0) + wa * wb
    patterns = {}
    for exact_intersection in overlap_binomial:
        dest = {}
        for subset_profile, coefficients in overlap_binomial.items():
            if any(a > b for a, b in zip(exact_intersection, subset_profile)):
                continue
            multiplier = (-1) ** (sum(subset_profile) - sum(exact_intersection))
            for a, b in zip(exact_intersection, subset_profile):
                multiplier *= comb(b, a)
            for total, weight in coefficients.items():
                dest[total] = dest.get(total, 0) + multiplier * weight
        assert all(w >= 0 for w in dest.values())
        for total, weight in dest.items():
            if weight:
                ones = tuple(a - 2 * b for a, b in zip(total, exact_intersection))
                assert min(ones) >= 0
                patterns[ones, exact_intersection] = patterns.get((ones, exact_intersection), 0) + weight
    assert sum(first_coef.values()) == 16 ** 5
    assert sum(patterns.values()) == 16 ** 10
    first, second = [F(0)] * 12, [F(0)] * 12
    for profile, weight in first_coef.items():
        first = [x + weight * y for x, y in zip(first, mixed(profile, (0, 0, 0)))]
    for (ones, twos), weight in patterns.items():
        second = [x + weight * y for x, y in zip(second, mixed(ones, twos))]
    assert first[0] == F(648, 1463) and second[0] == F(126404864, 474619145)
    assert first[1] == second[1] == 0
    return first, second, {'basis_group_coefficients': {str(k): v for k, v in sorted(first_coef.items())},
                           'pair_pattern_weights': {str(k): v for k, v in sorted(patterns.items())},
                           'pattern_count': len(patterns)}


def main():
    began = time.monotonic()
    ctx.prec = 256
    d1, d2, geometry = determinant_moments()
    bmoments = [F(0)] * 12
    # For the three-group Gram matrix, eigenvalues parallel/perpendicular to 1 are
    # 5a+9b+2c and 4b+12c, so 1'adj(A)1=5*(4b+12c)^4.
    for p in range(5):
        weight = 5 * comb(4, p) * 4 ** p * 12 ** (4 - p)
        bmoments = [x + weight * y for x, y in zip(bmoments, mixed((0, p, 4 - p), (0, 0, 0)))]
    assert bmoments[0] == F(216, 77)
    prior_path = ROOT.parent / 'astra_upper003' / 'verification.json'
    prior = json.loads(prior_path.read_text())
    A_lower = F(prior['constants']['A_lower'])
    B_upper = F(prior['constants']['B_upper'])
    c_lower = F(prior['constants']['c_lower'])
    c = 3 / (2 ** (arb(15) / 2) * arb.pi() ** (arb(5) / 2))
    ell = [A_lower * x - B_upper * y for x, y in zip(d1, d2)]
    candidates = []
    for eta in [F(1, 32), F(1, 16), F(1, 8)]:
        U = sum((F(comb(7, k)) * eta ** k * ell[k] for k in range(8)), F(0))
        V = sum((F(comb(11, k)) * eta ** k * bmoments[k] for k in range(12)), F(0))
        assert U > 0 and V > 0
        optimal_a = (ar(U) / ar(V)) ** (arb(1) / 4)
        tick = floor(float(optimal_a) * 1000)
        a = F(tick, 1000)
        assert ar(a) <= optimal_a < ar(a + F(1, 1000))
        value = F(11, 4) * U * a ** 7 - F(7, 4) * V * a ** 11
        ceiling = c * ar(U) ** (arb(11) / 4) / ar(V) ** (arb(7) / 4)
        candidates.append({'eta': str(eta), 'a_selected': str(a), 'optimal_a_ball': str(optimal_a),
                           'exact_corner_U': str(U), 'exact_corner_V': str(V),
                           'corner_polynomial_value': str(value), 'corner_g_value_Arb': str(c * ar(value)),
                           'all_a_corner_ceiling_Arb': str(ceiling),
                           'passes_necessary_corner_screen': bool(c * ar(value) > arb(1) / 27000),
                           'global_certificate_obtained': False})
    out = {'status': 'necessary_corner_checkpoint_only_not_a_global_bound',
           'new_g_bound': False, 'new_C5_upper_bound': False,
           'groups': {'weights': WEIGHTS, 'counts': COUNTS}, 'geometry': geometry,
           'moments': {'E_d_W_power': [str(x) for x in d1], 'E_d_squared_W_power': [str(x) for x in d2], 'E_b_corner_W_power': [str(x) for x in bmoments]},
           'constants': {'A_lower': str(A_lower), 'B_upper': str(B_upper), 'c_lower': str(c_lower)},
           'candidates': candidates,
           'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE, prior_path, ROOT / 'registration.md', Path(__file__)]},
           'runtime': {'python': sys.executable, 'arb_bits': ctx.prec, 'seconds': time.monotonic() - began}}
    (ROOT / 'corner_checkpoint.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'candidates': candidates, 'patterns': geometry['pattern_count'], 'runtime': out['runtime']}, indent=2))


if __name__ == '__main__':
    main()
