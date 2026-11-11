"""Local mpmath experiments for the r089 scalar Gaussian and explicit-source formulas.

These are numerical comparisons only. They do not certify physical norms or finite-N onset.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import json
import math
import multiprocessing as mp
import os
import time
from pathlib import Path

import mpmath as m

PRECISIONS = (80, 160)
K_VALUES = (16, 32, 64, 128, 256, 512, 1024)
ROOT_RATIOS = ('1', '1.25', '1.5', '1.75', '1.9')
QUANTILE_TYPES = ('1/N', '1/sqrt(NL)', 'K=.5sigma', 'K=.75sigma')
CLIPPING_MODES = ('P=.25', 'P=1', 'P=4', 'critical', 'quadratic_c=.1', 'quadratic_c=.2', 'quadratic_c=.3')


def _guard_digits(k):
    return math.ceil(k * math.log10(2) / 2) + 40


def _mp_nstr(x):
    return m.nstr(x, 65)


def _model_constants():
    mu = m.findroot(lambda u: u - 1 - m.exp(-u), (m.mpf('1.2'), m.mpf('1.4')))
    sigma = m.sqrt(mu)
    return mu, sigma


def _normal_parts(x):
    phi = m.exp(-x*x/2) / m.sqrt(2*m.pi)
    tail = m.erfc(x/m.sqrt(2)) / 2
    psi = phi - x*tail
    return phi, tail, psi


def _log_tail(x):
    return m.log(m.erfc(x/m.sqrt(2)) / 2)


def _solve_z(L, sigma):
    guess = m.sqrt(L)
    return m.findroot(lambda z: _log_tail(z) - m.log(sigma*z) + L/2,
                      (guess*m.mpf('0.8'), guess*m.mpf('1.1')))


def _format_root_row(k, ratio_s, wp):
    m.mp.dps = wp + _guard_digits(k)
    mu, sigma = _model_constants()
    L = m.mpf(k) * m.log(2)
    N = m.mpf(1 << k)
    sqrt_n = m.sqrt(N)
    ratio = m.mpf(ratio_s)
    z = _solve_z(L, sigma)
    w = ratio*m.sqrt(L) - z
    phi_z, _, psi_w = _normal_parts(z)[0], None, _normal_parts(w)[2]
    std_scale = sigma*psi_w/(sqrt_n*phi_z)
    log_D = m.log(sigma*sigma*psi_w/phi_z)

    def payoff(a):
        y = a*std_scale
        phi_zy, tail_zy, psi_zy = _normal_parts(z+y)
        _, tail_wy, psi_wy = _normal_parts(w-y)
        return -sqrt_n*tail_zy/sigma + (z+y) + psi_zy - psi_wy

    def derivative(a):
        y = a*std_scale
        _, tail_zy, _ = _normal_parts(z+y)
        _, tail_wy, _ = _normal_parts(w-y)
        phi_zy = m.exp(-(z+y)**2/2)/m.sqrt(2*m.pi)
        return std_scale*(sqrt_n*phi_zy/sigma + 1 - tail_zy - tail_wy)

    lo, hi = m.mpf(-2), m.mpf(2)
    flo, fhi = payoff(lo), payoff(hi)
    expansions = 0
    while flo*fhi > 0 and expansions < 60:
        if flo > 0:
            lo *= 2
            flo = payoff(lo)
        else:
            hi *= 2
            fhi = payoff(hi)
        expansions += 1
    if not flo <= 0 <= fhi:
        raise ArithmeticError(f'failed to bracket ramp root at k={k}, ratio={ratio_s}')
    initial_bracket = (lo, hi)

    # Twenty bisections provide a safe Newton starting point, not the final accuracy.
    for _ in range(20):
        mid = (lo+hi)/2
        fm = payoff(mid)
        if fm < 0:
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
    x = (lo+hi)/2
    target = m.power(10, -wp-8)
    iterations = 20
    for _ in range(2*wp+100):
        fx = payoff(x)
        dfx = derivative(x)
        if not dfx > 0:
            raise ArithmeticError(f'nonpositive ramp derivative at k={k}, ratio={ratio_s}')
        residual_step = abs(fx/dfx)
        if residual_step <= target:
            break
        candidate = x - fx/dfx
        if not lo < candidate < hi:
            candidate = (lo+hi)/2
        fc = payoff(candidate)
        if fc < 0:
            lo, flo = candidate, fc
        else:
            hi, fhi = candidate, fc
        x = candidate
        iterations += 1
    else:
        raise ArithmeticError(f'ramp root failed requested precision at k={k}, ratio={ratio_s}')
    residual_step = abs(payoff(x)/derivative(x))
    if residual_step > target:
        raise ArithmeticError(f'ramp residual/derivative too large at k={k}, ratio={ratio_s}')
    return {
        'working_precision': wp,
        'guarded_working_precision': m.mp.dps,
        'z': _mp_nstr(z), 'w': _mp_nstr(w), 'log_D': _mp_nstr(log_D),
        'normalized_root': _mp_nstr(x),
        'initial_bracket': [_mp_nstr(v) for v in initial_bracket],
        'final_bracket': [_mp_nstr(lo), _mp_nstr(hi)],
        'bracket_expansions': expansions, 'iterations': iterations,
        'residual_over_derivative': _mp_nstr(residual_step),
        'root_error_target': _mp_nstr(target),
    }


def root_job(job):
    k, ratio_s = job
    rows = [_format_root_row(k, ratio_s, wp) for wp in PRECISIONS]
    return {'k': k, 'N': str(1 << k), 'K_over_sigma': ratio_s, 'precision_rows': rows, 'pid': os.getpid()}


def _quantile_epsilon(kind, N, L, sigma):
    if kind == '1/N':
        return 1/N
    if kind == '1/sqrt(NL)':
        return 1/m.sqrt(N*L)
    rho = m.sqrt(L/N)
    if kind == 'K=.5sigma':
        return m.mpf('0.5')*sigma*rho
    if kind == 'K=.75sigma':
        return m.mpf('0.75')*sigma*rho
    raise ValueError(kind)


def _format_quantile_row(k, kind, wp):
    m.mp.dps = wp + _guard_digits(k)
    _, sigma = _model_constants()
    L = m.mpf(k)*m.log(2)
    N = m.mpf(1 << k)
    eps = _quantile_epsilon(kind, N, L, sigma)
    log_eps = m.log(eps)
    guess = m.sqrt(max(m.mpf(1), -2*log_eps))
    y = m.findroot(lambda x: _log_tail(x)-log_eps, (guess*m.mpf('.9'), guess*m.mpf('1.1')))
    z = _solve_z(L, sigma)
    shift = y-z
    loss = sigma*m.sqrt(N)*shift
    W = m.sqrt(N/L)
    normalized_loss = loss/W
    return {
        'working_precision': wp,
        'guarded_working_precision': m.mp.dps,
        'epsilon': _mp_nstr(eps),
        'log_epsilon': _mp_nstr(log_eps),
        'y': _mp_nstr(y),
        'z': _mp_nstr(z),
        'standardized_shift': _mp_nstr(shift),
        'gaussian_log_loss': _mp_nstr(loss),
        'gaussian_log_loss_over_W': _mp_nstr(normalized_loss),
    }


def quantile_job(job):
    k, kind = job
    rows = [_format_quantile_row(k, kind, wp) for wp in PRECISIONS]
    return {'k': k, 'N': str(1 << k), 'type': kind, 'precision_rows': rows, 'pid': os.getpid()}


def _e2(t, target):
    # E2(t) = sum_{j>=0} t^j/(j+2)!, evaluated by a stable recurrence.
    term = m.mpf(1)/2
    total = term
    j = 0
    while abs(term) > target*max(1, abs(total)):
        term *= t/(j+3)
        total += term
        j += 1
        if j > 10000:
            raise ArithmeticError('E2 series failed to converge')
    return total, j+1


def _solve_v(p, mu, target):
    r = m.exp(-mu)

    def residual(v):
        bp = mu+v/p
        t = -bp/p
        e2, _ = _e2(t, target/100)
        return v-p*r*m.expm1(-v/p)-bp*bp*e2

    left, right = mu/4, 3*mu/4
    fleft, fright = residual(left), residual(right)
    expansions = 0
    while fleft*fright > 0 and expansions < 40:
        width = right-left
        if fleft > 0:
            left -= width
            fleft = residual(left)
        else:
            right += width
            fright = residual(right)
        expansions += 1
    if not fleft <= 0 <= fright:
        raise ArithmeticError('failed to bracket v root')
    v = m.findroot(residual, (left, right), solver='anderson', tol=target, verify=True)
    err = abs(residual(v))
    if err > target*100:
        raise ArithmeticError('v root residual exceeds target')
    return v, err


def _clipping_parameters(mode, N, L, sigma):
    W = m.sqrt(N/L)
    rho = m.sqrt(L/N)
    if mode.startswith('P='):
        P = m.mpf(mode.split('=')[1])
        b = sigma/m.sqrt(2*P+1)
        p = b*W
        cap = sigma*(1+m.sqrt(2*P+1))*rho
        s = N**(-(P+1))
        normalized_by = W*m.log(L)
        return {'mode': mode, 'P': P, 'b': b, 'p': p, 'cap_excess': cap,
                's': s, 'normalization': normalized_by}
    if mode == 'critical':
        p = sigma*W
        cap = 2*sigma*rho
        s = 1/(N*L**m.mpf('.25'))
        return {'mode': mode, 'P': None, 'b': sigma, 'p': p, 'cap_excess': cap,
                's': s, 'normalization': W*m.log(L)}
    c = m.mpf(mode.split('=')[1])
    delta = N**(-m.mpf('.25'))
    b = sigma/m.sqrt(2*c)
    p = b/delta
    cap = delta
    s = m.exp(-c*N*delta*delta)
    return {'mode': mode, 'P': None, 'c': c, 'b': b, 'p': p,
            'cap_excess': cap, 'delta': delta, 's': s,
            'normalization': N*delta}


def _format_clipping_row(k, mode, wp):
    m.mp.dps = wp+_guard_digits(k)
    mu, sigma = _model_constants()
    L = m.mpf(k)*m.log(2)
    N = m.mpf(1 << k)
    params = _clipping_parameters(mode, N, L, sigma)
    p, b, cap, s = (params[q] for q in ('p', 'b', 'cap_excess', 's'))
    if p < 2:
        raise ArithmeticError(f'p<2 at k={k}, mode={mode}')
    if not (0 < s < 1):
        raise ArithmeticError(f's outside (0,1) at k={k}, mode={mode}')
    log_a = m.log1p(cap)+m.log1p(-s)
    if not (0 < log_a < m.log(2)):
        raise ArithmeticError(f'log(a) outside (0,log 2) at k={k}, mode={mode}')
    ln_phi = log_a+(N-1)*m.log1p(log_a)
    if not ln_phi > 0:
        raise ArithmeticError(f'ln Phi nonpositive at k={k}, mode={mode}')

    target = m.power(10, -wp-8)
    v, v_residual = _solve_v(p, mu, target)
    bp = mu+v/p
    _, sigma = _model_constants()
    z = _solve_z(L, sigma)
    C0 = 2*(1+m.log(2))
    log_s = m.log(s)
    high_gap = (N/p)*v + sigma*m.sqrt(N)*z + p*(m.log(C0)+m.log(p)-L-log_s)-ln_phi
    normalized_gap = high_gap/params['normalization']
    out = {
        'working_precision': wp,
        'guarded_working_precision': m.mp.dps,
        'p': _mp_nstr(p), 'b_parameter': _mp_nstr(b),
        'cap_excess': _mp_nstr(cap),
        's': None if log_s < -10000 else _mp_nstr(s),
        'threshold_value_omitted': bool(log_s < -10000), 'log_s': _mp_nstr(log_s),
        'log_a': _mp_nstr(log_a), 'ln_Phi_a': _mp_nstr(ln_phi),
        'C0': _mp_nstr(C0), 'v': _mp_nstr(v), 'b_p': _mp_nstr(bp),
        'v_equation_residual': _mp_nstr(v_residual),
        'v_root_error_target': _mp_nstr(target),
        'log_high_fee_minus_q0': _mp_nstr(high_gap),
        'normalized_gap': _mp_nstr(normalized_gap),
    }
    if 'P' in params:
        out['P'] = _mp_nstr(params['P']) if params['P'] is not None else None
    if 'c' in params:
        out['c'] = _mp_nstr(params['c'])
        out['delta'] = _mp_nstr(params['delta'])
    return out


def clipping_job(job):
    k, mode = job
    rows = [_format_clipping_row(k, mode, wp) for wp in PRECISIONS]
    return {'k': k, 'N': str(1 << k), 'mode': mode, 'precision_rows': rows, 'pid': os.getpid()}


def _check_finite_tree(obj, path='root'):
    if isinstance(obj, dict):
        for key, value in obj.items():
            _check_finite_tree(value, f'{path}.{key}')
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            _check_finite_tree(value, f'{path}[{i}]')
    elif isinstance(obj, str):
        if path.endswith(('N', 'mode', 'type', 'K_over_sigma')):
            return
        try:
            value = m.mpf(obj)
        except (ValueError, TypeError):
            return
        if not m.isfinite(value):
            raise ArithmeticError(f'nonfinite output at {path}: {obj}')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory')
    args = parser.parse_args()
    outdir = Path(args.output_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    root_jobs = [(k, ratio) for k in K_VALUES for ratio in ROOT_RATIOS]
    quantile_jobs = [(k, kind) for k in K_VALUES for kind in QUANTILE_TYPES]
    clipping_jobs = [(k, mode) for k in K_VALUES for mode in CLIPPING_MODES]
    with mp.get_context('spawn').Pool(12) as pool:
        root_rows = pool.map(root_job, root_jobs)
        quantile_rows = pool.map(quantile_job, quantile_jobs)
        clipping_rows = pool.map(clipping_job, clipping_jobs)
    _check_finite_tree(root_rows, 'root_rows')
    _check_finite_tree(quantile_rows, 'quantile_rows')
    _check_finite_tree(clipping_rows, 'clipping_rows')
    actual_workers = len({r['pid'] for r in root_rows+quantile_rows+clipping_rows})
    data = {
        'pass_all': True,
        'version': m.__version__,
        'configured_workers': 12,
        'actual_workers': actual_workers,
        'worker_pids': sorted({r['pid'] for r in root_rows+quantile_rows+clipping_rows}),
        'working_precision': list(PRECISIONS),
        'scope': 'Numerical scalar comparison of the Gaussian model and explicit high-source fee formulas only; no physical norm claim and no finite-N onset claim.',
        'counts': {'root_jobs': len(root_rows), 'quantile_jobs': len(quantile_rows),
                   'clipping_jobs': len(clipping_rows), 'total_jobs': len(root_rows)+len(quantile_rows)+len(clipping_rows)},
        'seconds': time.monotonic()-start,
        'root_rows': root_rows,
        'quantile_rows': quantile_rows,
        'clipping_rows': clipping_rows,
    }
    (outdir/'experiments89.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({k: v for k, v in data.items() if not k.endswith('_rows')}, indent=2))


if __name__ == '__main__':
    main()
