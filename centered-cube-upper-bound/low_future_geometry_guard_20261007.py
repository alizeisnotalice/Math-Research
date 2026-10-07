#!/usr/bin/env python3
"""Exact guards for new capture-count and low-Z column assertions only."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

OUT = Path(__file__).resolve().parent
PATH = OUT / 'low_future_geometry_guard_20261007_results.json'


def ratio(x):
    return f'{x.numerator}/{x.denominator}'


def digest_integer(x):
    raw = abs(x).to_bytes(max(1, (abs(x).bit_length() + 7) // 8), 'big')
    return dict(sign=1 if x > 0 else (-1 if x < 0 else 0),
                bit_length=abs(x).bit_length(), sha256_unsigned_big_endian=hashlib.sha256(raw).hexdigest())


def compact_fraction(x):
    return dict(numerator=digest_integer(x.numerator), denominator=digest_integer(x.denominator))


def log_core(z, terms):
    lower = 2 * sum((z ** (2 * k + 1) / (2 * k + 1) for k in range(terms)), F(0))
    upper = lower + 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    return lower, upper


def log_interval(x, terms=32):
    assert x > 0
    shift = 0
    while x < 1:
        x *= 2
        shift -= 1
    while x >= 2:
        x /= 2
        shift += 1
    lower, upper = log_core((x - 1) / (x + 1), terms)
    l2, u2 = log_core(F(1, 3), terms)
    if shift >= 0:
        return lower + shift * l2, upper + shift * u2
    return lower + shift * u2, upper + shift * l2


def main():
    assert not PATH.exists(), 'Frozen result refuses overwrite'
    rounds = []
    for n in [576, 1024, 4096]:
        alpha = isqrt(n)
        assert alpha * alpha == n
        ln_lo, ln_hi = log_interval(F(n + 2))
        ll_lo = log_interval(ln_lo)[0]
        ll_hi = log_interval(ln_hi)[1]
        ell_lo = 4 * F(alpha + 1, alpha) * ll_lo
        ell_hi = 4 * F(alpha + 1, alpha) * ll_hi
        records = []
        counts = dict(defect_large=0, not_large=0, ambiguous=0, baseline_large=0)
        for sigma in [F(2, n), F(1, 4), F(1, 2), F(3, 4), 1 - F(1, 2 * alpha)]:
            r = 1 - sigma
            log_sigma_lo, log_sigma_hi = log_interval(sigma)
            center_lo, center_hi = log_interval(sigma + r / F(3, 4))
            for phi in [F(3, 8), F(1, 2), F(3, 4)]:
                d_lo, d_hi = log_interval(sigma + r / phi)
                # Independent elementary log inequalities at the coordinate level.
                assert d_hi <= r * (1 / phi - 1) <= F(5, 3) * r
                assert log_sigma_hi <= -r
                for inside in [0, n // 4, 3 * n // 8, n // 2, 3 * n // 4, n]:
                    outside = n - inside
                    defect_lo = inside * d_lo + outside * log_sigma_lo
                    defect_hi = inside * d_hi + outside * log_sigma_hi
                    upper = r * (F(5, 3) * inside - outside)
                    assert defect_hi <= upper
                    baseline_lo = inside * center_lo + outside * log_sigma_lo
                    baseline_hi = inside * center_hi + outside * log_sigma_hi
                    assert baseline_lo <= defect_hi
                    if defect_lo > ell_hi:
                        classification = 'defect_large'
                        assert inside > F(3 * n, 8) + 3 * ell_lo / (8 * r)
                    elif defect_hi <= ell_lo:
                        classification = 'not_large'
                    else:
                        classification = 'ambiguous'
                    counts[classification] += 1
                    if baseline_lo > ell_hi:
                        counts['baseline_large'] += 1
                    records.append(dict(sigma=ratio(sigma), inside=inside, outside=outside,
                        phi_value_proxy=ratio(phi), defect_interval=[ratio(defect_lo), ratio(defect_hi)],
                        analytic_upper=ratio(upper), baseline_interval=[ratio(baseline_lo), ratio(baseline_hi)],
                        classification=classification))

        sigma = F(2, n)
        bound_a = F(3, 4) / (1 - sigma)
        assert bound_a <= F(4, 5)
        beta_bound = F(1, 2 ** alpha) + F(9, 10) ** n
        eps_lower = 1 / ln_hi ** 4
        assert beta_bound < eps_lower
        assert n * beta_bound < F(1, 16)
        survival = (1 - sigma) ** n
        assert survival >= F(1, 16)
        rounds.append(dict(n=n, alpha=alpha, capture_cases=len(records), classifications=counts,
            log_n_plus_two_interval=[ratio(ln_lo), ratio(ln_hi)],
            # Nested logarithm intervals have very large exact rational endpoints.
            # Hash their integers in binary; all comparisons above use full Fractions.
            defect_threshold_interval=[compact_fraction(ell_lo), compact_fraction(ell_hi)],
            records=records,
            low_Z_full_cube_guard=dict(inside_ratio_upper=ratio(bound_a),
                beta_average_upper=compact_fraction(beta_bound),
                epsilon_certified_lower=ratio(eps_lower),
                hard_component_survival=compact_fraction(survival),
                low_Z_less_than_epsilon=True, column_lower_coefficient='1/16',
                kernel_scope='Fixed-source free kernel only; no actual FIRST/nonconcentration/input instance.')))
    data = dict(status='PASS_THREE_ROUNDS_EXACT_NEW_LOW_FUTURE_GEOMETRY_GUARDS',
        rounds=rounds, capture_cases_total=sum(x['capture_cases'] for x in rounds),
        arithmetic='Fraction and positive log series with explicit rational remainder; no floating pass predicate.',
        ratio_proxy_scope='Coordinate phi proxies within [3/8,3/4], not numerical samples of the original phi.',
        scope='Only new capture-count consequence and entire-initial-cube low-Z coefficient bounds; no repeated Jensen tests, no actual FIRST counterexample or source fee.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with PATH.open('x') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(dict(status=data['status'], capture_cases=data['capture_cases_total'],
        rounds=[dict(n=x['n'], counts=x['classifications']) for x in rounds])))


if __name__ == '__main__':
    main()
