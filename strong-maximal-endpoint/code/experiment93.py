"""r093 exact rational checks for selected hyperbolic-halo inequalities.

The first batch certifies s_p > N b_p on a finite (N,p) grid using the
Gamma(N,1) integral bound. The second checks the stated finite polynomial
family. Neither batch proves a general H_N bound or computes sharp s_p values.
"""
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

import mpmath as m

HERE = Path(__file__).resolve().parent
T_FALLBACK = Path('/Volumes/data/CodexProjects/T')
PROJECT_CODE = T_FALLBACK/'Math-Research-strong-maximal/strong-maximal-endpoint/code'
CODE = next((path for path in (HERE/'code', HERE, PROJECT_CODE)
             if (path/'sliver79.py').is_file()), None)
if CODE is None:
    raise ImportError('could not locate sliver79.py')
sys.path.insert(0, str(CODE))
from sliver79 import logi

N_VALUES = (2, 3, 4, 8, 16, 32, 64)
P_VALUES = (Q(9, 10), Q(19, 20), Q(49, 50))
A_VALUES = (Q(0), Q(1, 4), Q(1), Q(4))
PRECISIONS = (80, 160)
WORKERS = 12
SIG_DIGITS = 40


def fraction_digest(value):
    assert isinstance(value, Q)
    canonical = f'{value.numerator}/{value.denominator}'.encode('ascii')
    return {
        'sha256': hashlib.sha256(canonical).hexdigest(),
        'numerator_bits': abs(value.numerator).bit_length(),
        'denominator_bits': value.denominator.bit_length(),
    }


def _pow10(exponent):
    return Q(10**exponent) if exponent >= 0 else Q(1, 10**(-exponent))


def floor_log10_fraction(value):
    """Exact floor(log10(value)) for a positive Fraction."""
    assert value > 0
    # A rational bit-length estimate followed by exact comparisons avoids any
    # floating-point dependence, even for tiny values.
    e = ((value.numerator.bit_length()-value.denominator.bit_length())*30103)//100000
    while value < _pow10(e):
        e -= 1
    while value >= _pow10(e+1):
        e += 1
    return e


def scientific_enclosure(value, digits=SIG_DIGITS):
    """Compact exact decimal enclosure: [lo,hi]*10**(e-digits+1)."""
    assert value >= 0
    if value == 0:
        return {
            'lower_significand_integer': '0',
            'upper_significand_integer': '0',
            'decimal_exponent': 0,
            'significand_scale_power': digits-1,
            'exact_zero': True,
        }
    e = floor_log10_fraction(value)
    scale_power = e-digits+1
    scaled = value/_pow10(scale_power)
    lo = scaled.numerator//scaled.denominator
    hi = lo if scaled.denominator == 1 else lo+1
    assert Q(lo, 10**(digits-1))*_pow10(e) <= value
    assert value <= Q(hi, 10**(digits-1))*_pow10(e)
    return {
        'lower_significand_integer': str(lo),
        'upper_significand_integer': str(hi),
        'decimal_exponent': e,
        'significand_scale_power': digits-1,
        'exact_zero': False,
    }


def log_interval(value):
    lo, hi = logi(value)
    return [str(lo), str(hi)]


def _mpf_fraction(value):
    return m.mpf(value.numerator)/value.denominator


def batch1_job(arg):
    n, p_string = arg
    p = Q(p_string)
    b_argument = Q(2)/p-1
    b_lo, b_hi = logi(b_argument)
    assert b_lo <= b_hi and b_argument > 0

    # Gamma(N,1) CDF at Nb is at most (Nb)^N/N! because exp(-t)<=1.
    cdf_upper = (n*b_hi)**n/Q(math.factorial(n))
    threshold = 1-p
    margin = threshold-cdf_upper
    assert margin > 0, (n, p_string, cdf_upper, threshold)
    assert cdf_upper < threshold

    numeric = []
    for wp in PRECISIONS:
        with m.workdps(wp):
            b_num = m.log(2/_mpf_fraction(p)-1)
            log_cdf = n*m.log(n*b_num)-m.loggamma(n+1)
            cdf_num = m.exp(log_cdf)
            numeric.append({
                'working_precision': wp,
                'b': m.nstr(b_num, 65),
                'log_gamma_cdf_upper_at_exact_b': m.nstr(log_cdf, 65),
                'gamma_cdf_upper_at_exact_b': m.nstr(cdf_num, 65),
                'diagnostic_only': True,
            })

    return {
        'N': n,
        'p': p_string,
        'b_argument_2_over_p_minus_1': str(b_argument),
        'b_interval': [str(b_lo), str(b_hi)],
        'gamma_cdf_upper_formula': '(N*b_hi)^N/N!',
        'gamma_cdf_upper_log_interval': log_interval(cdf_upper),
        'gamma_cdf_upper_exact_fraction_sha256': fraction_digest(cdf_upper),
        'gamma_cdf_upper_tight_decimal_enclosure': scientific_enclosure(cdf_upper),
        'threshold_1_minus_p': str(threshold),
        'strict_margin_exact_fraction_sha256': fraction_digest(margin),
        'strict_margin_lower_decimal_bound': {
            key: value for key, value in scientific_enclosure(margin).items()
            if key != 'upper_significand_integer'
        },
        'strict_margin_log_interval': log_interval(margin),
        'exact_fraction_comparison': 'gamma_cdf_upper < 1-p',
        'strict_margin_positive_exact': True,
        'integral_bound_uses_exp_minus_t_le_1': True,
        'precision_rows': numeric,
        'pid': os.getpid(),
    }


def polynomial(n, t):
    assert t >= 0
    total = Q(0)
    term = Q(1)
    total += term
    for k in range(1, n):
        term *= t/Q(k)
        total += term
    return total


def batch2_job(n):
    rows = []
    for a in A_VALUES:
        for t in (Q(0), a/2, a, a+1, a+4):
            pt = polynomial(n, t)
            if t >= a:
                numerator = polynomial(n, t-a)
                ratio = numerator/pt
                assert 0 <= ratio <= 1
                if a == 0:
                    assert ratio == 1
                else:
                    assert ratio < 1
                rows.append({
                    'N': n,
                    'a': str(a),
                    'T': str(t),
                    'case': 'T_ge_a_polynomial_ratio',
                    'exact_ray_formula': 'P_{N-1}(T-a)/P_{N-1}(T)',
                    'exact_ray_upper_fraction_sha256': fraction_digest(ratio),
                    'exact_ray_upper_tight_decimal_enclosure': scientific_enclosure(ratio),
                    'exact_ray_le_one': True,
                    'exact_ray_strictly_less_than_one': ratio < 1,
                    'exact_fraction_le_one': ratio <= 1,
                    'equality_at_one': ratio == 1,
                    'outside_halo_upper_fraction_sha256': None,
                    'outside_halo_upper_tight_decimal_enclosure': None,
                    'analytic_exp_factor_upper': None,
                    'analytic_exp_factor_is_equality': None,
                })
            else:
                reciprocal = 1/pt
                assert pt >= 1 and reciprocal <= 1
                # Exact outside-halo factor is exp(T-a); since T<a it is
                # strictly less than 1. Store the rational upper envelope
                # 1/P(T), never mistake it for the exponential ratio itself.
                rows.append({
                    'N': n,
                    'a': str(a),
                    'T': str(t),
                    'case': 'T_lt_a_exponential_times_reciprocal',
                    'exact_ray_formula': 'exp(T-a)/P_{N-1}(T)',
                    'exact_ray_upper_fraction_sha256': None,
                    'exact_ray_upper_tight_decimal_enclosure': None,
                    'exact_ray_le_one': True,
                    'exact_ray_strictly_less_than_one': True,
                    'equality_at_one': False,
                    'outside_halo_upper_fraction_sha256': fraction_digest(reciprocal),
                    'outside_halo_upper_tight_decimal_enclosure': scientific_enclosure(reciprocal),
                    'outside_halo_upper_formula': '1/P_{N-1}(T)',
                    'analytic_exp_factor_upper': '1',
                    'analytic_exp_factor_is_equality': False,
                    'analytic_reason': 'T-a<0 implies exp(T-a)<1',
                })
    assert len(rows) == len(A_VALUES)*5
    return {'N': n, 'cases': rows, 'pid': os.getpid()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory')
    args = parser.parse_args()
    outdir = Path(args.output_directory)
    outdir.mkdir(parents=True, exist_ok=True)

    start = time.monotonic()
    jobs1 = [(n, str(p)) for n in N_VALUES for p in P_VALUES]
    with mp.get_context('spawn').Pool(WORKERS) as pool:
        batch1 = pool.map(batch1_job, jobs1)
        batch2 = pool.map(batch2_job, N_VALUES)

    # Strip process IDs from certificate rows; preserve distinct worker IDs as
    # execution metadata only.
    worker_ids = {row['pid'] for row in batch1+batch2}
    for row in batch1:
        del row['pid']
    for row in batch2:
        del row['pid']
    assert len(batch1) == 21 and len(batch2) == 7
    all_polynomial = [case for group in batch2 for case in group['cases']]
    assert len(all_polynomial) == 140
    assert all(row['strict_margin_positive_exact'] for row in batch1)
    assert all(row.get('exact_fraction_le_one', row.get('exact_ray_le_one'))
               for row in all_polynomial)
    assert all(row['exact_ray_le_one'] for row in all_polynomial)
    assert all(row['exact_ray_strictly_less_than_one'] == (row['a'] != '0')
               for row in all_polynomial)
    assert all(row['equality_at_one'] == (row['a'] == '0')
               for row in all_polynomial)
    assert all(row['case'] != 'T_lt_a_exponential_times_reciprocal' or
               row['analytic_exp_factor_is_equality'] is False
               for row in all_polynomial)
    assert all(row['equality_at_one'] for row in all_polynomial
               if row['a'] == '0')

    result = {
        'pass_all': True,
        'configured_workers': WORKERS,
        'actual_workers_observed_in_polynomial_batch': len(worker_ids),
        'workers_seen': sorted(worker_ids),
        'N_values': list(N_VALUES),
        'p_values': [str(p) for p in P_VALUES],
        'a_values': [str(a) for a in A_VALUES],
        'precision_rows': list(PRECISIONS),
        'batch1_gamma_cdf_cases': len(batch1),
        'batch2_polynomial_cases': len(all_polynomial),
        'batch1': batch1,
        'batch2': batch2,
        'seconds': time.monotonic()-start,
        'scope': ('Exact finite rational certificates for s_p>N*b_p on the listed grid '
                  'and the stated finite polynomial family only; no general H_N bound, '
                  'no general hyperbolic maximal theorem, and no sharp numerical s_p claim.'),
        'general_H_N_upper_bound_claimed': False,
        'sharp_s_p_value_claimed': False,
    }
    (outdir/'certificates93.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('batch1', 'batch2')}, indent=2))


if __name__ == '__main__':
    main()
