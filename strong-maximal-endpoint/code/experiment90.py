"""Exact/numeric r090 quadratic-price normal-form scale experiment.

The rational intervals certify the scalar Taylor remainder bounds. The mpmath
rows are numerical diagnostics and do not establish the physical witness theorem.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import json
import math
import multiprocessing as mp
import time
from fractions import Fraction as Q
from pathlib import Path
import sys

import mpmath as m

HERE = Path(__file__).resolve().parent
WORK_PUBLIC_FALLBACK = Path('/Volumes/data/CodexProjects/T/Math-Research-strong-maximal/strong-maximal-endpoint/code')
PUBLIC_CODE = HERE if (HERE/'verify85.py').is_file() else WORK_PUBLIC_FALLBACK
sys.path.insert(0, str(PUBLIC_CODE))
from verify85 import constants85, sqrt_interval
from sliver79 import logi
from experiment88 import log_large_fraction

K_VALUES = (16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024)
PRECISIONS = (80, 160)


def _mul_interval(a, b):
    products = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return min(products), max(products)


def _ceil_fraction(q):
    return -((-q.numerator)//q.denominator)


def _sqrt_log_interval(log_interval):
    assert log_interval[0] > 0
    return sqrt_interval(log_interval)


def _exact_scales(k, mu_iv):
    n = 1 << k
    log2_iv = logi(Q(2))
    L_iv = (k*log2_iv[0], k*log2_iv[1])
    muL_iv = (mu_iv[0]*L_iv[0], mu_iv[1]*L_iv[1])
    sqrt_muL = sqrt_interval(muL_iv)
    sqrt_n = 1 << (k//2)  # all listed k are even
    x_iv = (2*sqrt_muL[0]/sqrt_n, 2*sqrt_muL[1]/sqrt_n)
    assert 0 < x_iv[0] <= x_iv[1] < Q(1, 4)

    # D3 decreases on the required interval. Use the upper endpoint in both
    # decreasing factors to obtain a directed lower bound.
    u_hi = logi(1+x_iv[1])[1]
    h_hi = (2*u_hi*u_hi+7*u_hi+7)/(1+u_hi)**3
    d3_lower = (2+(n-1)*h_hi)/(1+x_iv[1])**3
    d3_upper = Q(7*n-5)
    rem_lo = x_iv[0]**3*d3_lower/6
    rem_hi = x_iv[1]**3*d3_upper/6
    assert 0 < rem_lo <= rem_hi
    rem_log_lo = log_large_fraction(rem_lo)[0]
    rem_log_hi = log_large_fraction(rem_hi)[1]

    cubic_lo = x_iv[0]**3*d3_upper/6
    cubic_hi = x_iv[1]**3*d3_upper/6
    cubic_log_iv = (log_large_fraction(cubic_lo)[0], log_large_fraction(cubic_hi)[1])

    # h=L N^(-1/4); k is divisible by four for every requested dimension.
    assert k % 4 == 0
    n_quarter = 1 << (k//4)
    h_iv = (L_iv[0]/n_quarter, L_iv[1]/n_quarter)
    sigma_iv = sqrt_interval(mu_iv)
    sqrtL_iv = _sqrt_log_interval(L_iv)
    grid_ratio_upper = 10*sigma_iv[1]*sqrtL_iv[1]*n_quarter/L_iv[0]
    grid_count_upper = _ceil_fraction(grid_ratio_upper)+5
    assert grid_count_upper < n
    return {
        'N': n, 'L_interval': L_iv, 'x_max_interval': x_iv,
        'remainder_log_interval': (rem_log_lo, rem_log_hi),
        'cubic_bound_log_interval': cubic_log_iv,
        'h_interval': h_iv, 'grid_count_upper': grid_count_upper,
    }


def _mp_root_mu(mu_iv):
    # The exact interval brackets the unique root; high precision findroot is
    # used only to evaluate the independent scalar integral.
    start = (m.mpf(mu_iv[0].numerator)/mu_iv[0].denominator +
             m.mpf(mu_iv[1].numerator)/mu_iv[1].denominator)/2
    return m.findroot(lambda u: u-1-m.exp(-u), start)


def _log_normal_tail(z):
    return m.log(m.erfc(z/m.sqrt(2))/2)


def _numeric_row(k, wp, mu_iv, exact):
    guard = math.ceil(k*math.log10(2)/2)+40
    m.mp.dps = wp+guard
    N = m.mpf(1 << k)
    L = m.mpf(k)*m.log(2)
    mu = _mp_root_mu(mu_iv)
    sigma = m.sqrt(mu)
    x = 2*m.sqrt(mu*L/N)
    if not (x > 0 and x < m.mpf(1)/4):
        raise ArithmeticError(f'x outside Taylor domain for k={k}')

    def third_derivative(t):
        u = m.log1p(x*t)
        hval = (2*u*u+7*u+7)/(1+u)**3
        return (2+(N-1)*hval)/m.exp(3*u)

    integral = m.quad(lambda t: (1-t)**2*third_derivative(t), [0, m.mpf('0.5'), 1])
    remainder = x**3*integral/2
    cubic = (7*N-5)*x**3/6
    if not (remainder > 0 and cubic > 0 and remainder <= cubic*(1+m.mpf('1e-40'))):
        raise ArithmeticError(f'quadratic remainder/cubic bound check failed at k={k}')
    log_remainder = m.log(remainder)
    log_cubic = m.log(cubic)
    ratio = remainder/cubic

    # Asymptotic scale markers are scalar diagnostics only.
    z0 = m.sqrt(L)
    z = m.findroot(lambda zz: _log_normal_tail(zz)-m.log(sigma*zz)+L/2,
                   (z0*m.mpf('.8'), z0*m.mpf('1.1')))
    h = L/N**m.mpf('.25')
    delta0 = sigma*z/m.sqrt(N)
    delta0_over_h = delta0/h
    normalized_identity = (L**2/(h*m.sqrt(N)))/h
    if abs(normalized_identity-1) > m.power(10, -wp-10):
        raise ArithmeticError(f'normalized scale identity failed at k={k}')

    # The directed logarithmic interval should enclose the independently
    # integrated remainder; it is deliberately not expected to be narrow.
    log_lo = m.mpf(exact['remainder_log_interval'][0].numerator)/exact['remainder_log_interval'][0].denominator
    log_hi = m.mpf(exact['remainder_log_interval'][1].numerator)/exact['remainder_log_interval'][1].denominator
    if not log_lo <= log_remainder <= log_hi:
        raise ArithmeticError(f'numerical integral lies outside exact log remainder bounds at k={k}, wp={wp}')

    return {
        'working_precision': wp,
        'guarded_working_precision': m.mp.dps,
        'mu': m.nstr(mu, 65), 'sigma': m.nstr(sigma, 65),
        'L': m.nstr(L, 65), 'x_max': m.nstr(x, 65),
        'log_remainder': m.nstr(log_remainder, 65),
        'log_cubic_bound': m.nstr(log_cubic, 65),
        'remainder_to_cubic_bound': m.nstr(ratio, 65),
        'z': m.nstr(z, 65), 'h': m.nstr(h, 65),
        'delta0': m.nstr(delta0, 65),
        'delta0_over_h': m.nstr(delta0_over_h, 65),
        'normalized_identity_value': m.nstr(normalized_identity, 65),
        'grid_count_upper': exact['grid_count_upper'],
    }


def job(arg):
    k, mu_strings = arg
    mu_iv = tuple(Q(s) for s in mu_strings)
    exact = _exact_scales(k, mu_iv)
    rows = [_numeric_row(k, wp, mu_iv, exact) for wp in PRECISIONS]
    return {
        'k': k, 'N': str(exact['N']),
        'x_max_interval': [str(x) for x in exact['x_max_interval']],
        'remainder_log_interval': [str(x) for x in exact['remainder_log_interval']],
        'uniform_remainder_log_upper': str(exact['remainder_log_interval'][1]),
        'uniform_remainder_is_error_bound': True,
        'cubic_bound_log_interval': [str(x) for x in exact['cubic_bound_log_interval']],
        'h_interval': [str(x) for x in exact['h_interval']],
        'normalized_scale_identity_exact': '1',
        'grid_count_upper': exact['grid_count_upper'],
        'grid_enumerated': False,
        'precision_rows': rows,
        'pid': os.getpid(),
    }


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory')
    args = parser.parse_args()
    outdir = Path(args.output_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    constants = constants85()
    mu_iv = tuple(Q(s) for s in constants['mu_interval'])
    if not (len(mu_iv) == 2 and mu_iv[0] < mu_iv[1]):
        raise ArithmeticError('invalid certified mu interval')
    start = time.monotonic()
    with mp.get_context('spawn').Pool(12) as pool:
        rows = pool.map(job, [(k, constants['mu_interval']) for k in K_VALUES])
    actual_workers = len({r['pid'] for r in rows})
    data = {
        'pass_all': True,
        'version': m.__version__,
        'configured_workers': 12,
        'actual_workers': actual_workers,
        'worker_pids': sorted({r['pid'] for r in rows}),
        'working_precision': list(PRECISIONS),
        'constants': {'mu_interval': constants['mu_interval'],
                      'sigma_interval': [str(x) for x in sqrt_interval(mu_iv)]},
        'scope': 'Exact scalar Taylor-remainder bounds and numerical scale diagram only; no numerical proof of the physical witness theorem or finite-N physical onset.',
        'counts': {'dimensions': len(rows), 'numeric_precision_rows': sum(len(r['precision_rows']) for r in rows),
                   'grid_points_enumerated': 0},
        'seconds': time.monotonic()-start,
        'rows': rows,
    }
    if len(rows) != len(K_VALUES) or any(r['grid_enumerated'] for r in rows):
        raise ArithmeticError('dimension or non-enumeration check failed')
    (outdir/'experiments90.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({k: v for k, v in data.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
