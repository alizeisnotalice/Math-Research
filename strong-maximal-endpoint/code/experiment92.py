"""r092 scalar rank-shift and positive-part diagnostics.

Only deterministic scalar calculations are performed. The output is not a
certificate for the unknown r021 constant and is not a physical onset claim.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import json
import math
import multiprocessing as mp
import time
from pathlib import Path
import sys

import mpmath as m

HERE = Path(__file__).resolve().parent
T_FALLBACK = Path('/Volumes/data/CodexProjects/T')
CANDIDATES = (HERE/'code', HERE,
              T_FALLBACK/'Math-Research-strong-maximal/strong-maximal-endpoint/code')
CODE = next((p for p in CANDIDATES if (p/'verify85.py').is_file()), None)
if CODE is None:
    raise ImportError('could not locate verify85.py in local code or T fallback')
sys.path.insert(0, str(CODE))
from verify85 import constants85

K_VALUES = (16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
B_VALUES = ('1/4', '1', '4', '16')
D_VALUES = (0, 1, 4)
PRECISIONS = (80, 160)


def _mp_mu(mu_interval):
    endpoints = [m.mpf(s.split('/')[0])/m.mpf(s.split('/')[1]) if '/' in s else m.mpf(s)
                 for s in mu_interval]
    start = (endpoints[0]+endpoints[1])/2
    return m.findroot(lambda u: u-1-m.exp(-u), start)


def _normal_parts(x):
    phi = m.exp(-x*x/2)/m.sqrt(2*m.pi)
    tail = m.erfc(x/m.sqrt(2))/2
    cdf = 1-tail
    return phi, tail, cdf


def _one_precision(k, b_string, wp, mu_interval):
    guard = math.ceil(k*math.log10(2)/2)+80
    m.mp.dps = wp+guard
    N_integer = 1 << k
    N = m.mpf(N_integer)
    L = m.mpf(k)*m.log(2)
    mu = _mp_mu(mu_interval)
    sigma = m.sqrt(mu)
    sqrtN = m.sqrt(N)
    B = m.mpf(b_string)

    def z_residual(z):
        return (m.log(m.erfc(z/m.sqrt(2))/2)-m.log(sigma*z)+L/2)

    z_guess = m.sqrt(L)
    z = m.findroot(z_residual, (z_guess*m.mpf('.8'), z_guess*m.mpf('1.1')))
    phi_z, tail_z, _ = _normal_parts(z)
    z_derivative = -phi_z/tail_z-1/z
    z_log_resid = abs(z_residual(z))
    z_normalized_resid = abs(z_residual(z)/z_derivative)

    tau = sigma*z/sqrtN
    target_tail = tau+B*L/N
    if not (0 < target_tail < m.mpf('0.5')):
        raise ArithmeticError(f'B-tail target outside (0,1/2) for k={k}, B={b_string}')

    def zb_residual(zb):
        return m.log(m.erfc(zb/m.sqrt(2))/2)-m.log(target_tail)

    # The B-tail target is above tau, so z_B is just below z; this local
    # two-point solve retains the rank signal instead of subtracting mu*N ranks.
    zB = m.findroot(zb_residual, (z-1, z))
    phi_zB, tail_zB, _ = _normal_parts(zB)
    zB_derivative = -phi_zB/tail_zB
    zB_log_resid = abs(zb_residual(zB))
    zB_normalized_resid = abs(zb_residual(zB)/zB_derivative)
    target_resid = m.power(10, -wp)
    if max(z_log_resid, zB_log_resid) >= target_resid:
        raise ArithmeticError(f'root log residual exceeds target for k={k}, B={b_string}, wp={wp}')

    # Form the local rank difference before scaling. Never subtract quantities
    # of size mu*N; guarded precision resolves the close normal quantiles.
    rank_shift = sigma*sqrtN*(z-zB)
    rank_shift_over_B = rank_shift/B
    if rank_shift <= 0:
        raise ArithmeticError(f'nonpositive rank shift for k={k}, B={b_string}, wp={wp}')

    linear_rows = []
    for Dvalue in D_VALUES:
        a = z+Dvalue*L/(sigma*sqrtN)
        phi_a, tail_a, cdf_a = _normal_parts(a)
        psi_a = phi_a-a*tail_a
        correction = sigma*sqrtN*psi_a
        normalized_excess = Dvalue+correction/L
        second_moment_scale = mu*((a*a+1)*cdf_a+a*phi_a)/L
        linear_rows.append({
            'D': Dvalue,
            'a': m.nstr(a, 65),
            'phi_a': m.nstr(phi_a, 65),
            'barPhi_a': m.nstr(tail_a, 65),
            'psi_a': m.nstr(psi_a, 65),
            'stoploss_correction': m.nstr(correction, 65),
            'N_times_normalized_excess': m.nstr(Dvalue*L+correction, 65),
            'normalized_excess_error_over_L': m.nstr(normalized_excess, 65),
            'second_moment_over_L_scale': m.nstr(second_moment_scale, 65),
        })

    return {
        'working_precision': wp,
        'guarded_working_precision': m.mp.dps,
        'N': str(N_integer), 'L': m.nstr(L, 65),
        'mu': m.nstr(mu, 65), 'sigma': m.nstr(sigma, 65),
        'z': m.nstr(z, 65), 'zB': m.nstr(zB, 65),
        'tau': m.nstr(tau, 65), 'B': m.nstr(B, 65),
        'target_tail_tau_plus_BL_over_N': m.nstr(target_tail, 65),
        'z_log_residual': m.nstr(z_log_resid, 65),
        'z_residual_over_derivative': m.nstr(z_normalized_resid, 65),
        'zB_log_residual': m.nstr(zB_log_resid, 65),
        'zB_residual_over_derivative': m.nstr(zB_normalized_resid, 65),
        'root_log_residual_target': m.nstr(target_resid, 65),
        'rank_shift': m.nstr(rank_shift, 65),
        'rank_shift_over_B': m.nstr(rank_shift_over_B, 65),
        'linear_positive_part_rows': linear_rows,
    }


def job(arg):
    k, B, mu_interval = arg
    rows = [_one_precision(k, B, wp, mu_interval) for wp in PRECISIONS]
    return {'k': k, 'N': str(1 << k), 'B': B,
            'precision_rows': rows, 'pid': os.getpid()}


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory')
    args = parser.parse_args()
    outdir = Path(args.output_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    constants = constants85()
    mu_interval = constants['mu_interval']
    jobs = [(k, B, mu_interval) for k in K_VALUES for B in B_VALUES]
    start = time.monotonic()
    with mp.get_context('spawn').Pool(12) as pool:
        rows = pool.map(job, jobs)
    actual_workers = len({r['pid'] for r in rows})
    data = {
        'pass_all': True,
        'version': m.__version__,
        'configured_workers': 12,
        'actual_workers': actual_workers,
        'worker_pids': sorted({r['pid'] for r in rows}),
        'working_precision': list(PRECISIONS),
        'precision_guard': 'wp + ceil(log10(N)/2) + 80',
        'constants': {'mu_interval': mu_interval},
        'scope': 'Deterministic scalar rank-shift and positive-part diagnostics only; B is not an authentication of the unknown r021 constant, the second-moment row is not an inverse-Phi expectation correction, and no finite-N physical onset is claimed.',
        'counts': {'k_values': len(K_VALUES), 'B_values': len(B_VALUES),
                   'root_jobs': len(rows), 'precision_rows': len(rows)*len(PRECISIONS),
                   'linear_positive_part_cases': len(rows)*len(PRECISIONS)*len(D_VALUES),
                   'high_dimensional_arrays_allocated': 0},
        'seconds': time.monotonic()-start,
        'rows': rows,
    }
    if len(rows) != 36:
        raise ArithmeticError('expected 36 scalar jobs')
    (outdir/'experiments92.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({k: v for k, v in data.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
