#!/usr/bin/env python3
"""Finite CDF/Bellman comparison for predictable-prefix feedback.

All probabilities and state values are numerical float64 diagnostics. The
fixed-action terminal law is evaluated as an exact-form Gamma--Binomial sum;
the Bellman step uses analytic exponential convolution of the piecewise-linear
state interpolant and positive linear interpolation for action shifts.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import time
from pathlib import Path

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import expit, gammainc
from scipy.stats import binom, norm

ROOT = Path(__file__).resolve().parent
NS = (32, 64, 128, 256, 512)
DXS = (1 / 16, 1 / 32, 1 / 64)
FRACTIONS = (1 / 4, 1 / 2, 3 / 4)
WORKERS = 12
ACTION_GRID = np.arange(0.0, 4.0 + 1e-12, 1 / 16)


def constants():
    mu = brentq(lambda x: x - 1.0 - math.exp(-x), 1.0, 2.0,
                xtol=2e-15, rtol=4 * np.finfo(float).eps)
    return mu, math.sqrt(mu)


MU, SIGMA = constants()


def normal_tail_root(n: int) -> tuple[float, float]:
    root = brentq(lambda z: norm.sf(z) - SIGMA * z / math.sqrt(n),
                  0.0, 8.0, xtol=1e-14, rtol=4 * np.finfo(float).eps)
    return root, SIGMA * root / math.sqrt(n)


def action_c(n: int) -> tuple[float, float, float]:
    z, alpha = normal_tail_root(n)
    c = MU - SIGMA * z / (2.0 * math.sqrt(n))
    if not (0.0 <= c <= 4.0):
        raise ArithmeticError(f'c_N outside action range for N={n}: {c}')
    return c, z, alpha


def gamma_binomial_cdf_vector(x: np.ndarray, n: int, c: float) -> np.ndarray:
    """CDF of Gamma(n,1)+c*Binomial(n, logistic(-c))."""
    q = float(expit(-c))
    j = np.arange(n + 1)
    weights = np.asarray(binom.pmf(j, n, q), dtype=np.float64)
    out = np.zeros_like(x, dtype=np.float64)
    for jj, weight in enumerate(weights):
        if weight:
            out += weight * gammainc(float(n), np.maximum(x - c * jj, 0.0))
    out[0] = 0.0
    return out


def gamma_binomial_cdf_scalar(x: float, n: int, c: float) -> float:
    q = float(expit(-c))
    j = np.arange(n + 1)
    weights = binom.pmf(j, n, q)
    vals = gammainc(float(n), np.maximum(x - c * j, 0.0))
    return float(np.dot(weights, vals))


def exp_convolution_piecewise_linear(f: np.ndarray, dx: float) -> np.ndarray:
    """G_i=integral_0^{x_i} exp(-(x_i-y)) F_linear(y) dy exactly per cell."""
    decay = math.exp(-dx)
    area = -math.expm1(-dx)
    first_moment = 1.0 - (1.0 + dx) * decay
    forcing = np.zeros_like(f)
    slope = np.diff(f) / dx
    forcing[1:] = area * f[1:] - first_moment * slope
    conv = np.empty_like(f)
    conv[0] = 0.0
    for i in range(1, len(f)):
        conv[i] = decay * conv[i - 1] + forcing[i]
    conv[0] = 0.0
    return conv


def shifted_grid_indices(size: int, dx: float, controls: np.ndarray):
    ix = np.arange(size, dtype=np.int64)[None, :]
    offsets = controls / dx
    nearest = np.rint(offsets)
    aligned = np.abs(offsets - nearest) < 2e-12
    whole = np.floor(offsets).astype(np.int64)
    whole[aligned] = nearest[aligned].astype(np.int64)
    frac = offsets - whole
    hi = ix - whole[:, None]
    lo = hi - 1
    valid_hi = hi >= 0
    valid_lo = lo >= 0
    return np.clip(hi, 0, size - 1), np.clip(lo, 0, size - 1), valid_hi, valid_lo, frac


def interp_quantile(x: np.ndarray, f: np.ndarray, probability: float) -> tuple[float, int, float]:
    if not (f[0] <= probability <= f[-1]):
        raise ArithmeticError(f'quantile target outside grid CDF: target={probability}, Fends=({f[0]},{f[-1]})')
    idx = int(np.searchsorted(f, probability, side='left'))
    if idx <= 0:
        return float(x[0]), 0, 0.0
    if idx >= len(f):
        return float(x[-1]), len(f) - 1, 0.0
    den = float(f[idx] - f[idx - 1])
    if den <= 0:
        raise ArithmeticError('flat CDF segment at target quantile')
    frac = (probability - float(f[idx - 1])) / den
    return float(x[idx - 1] + frac * (x[idx] - x[idx - 1])), idx, frac


def build_cases(pilot: bool):
    if pilot:
        return ([(32, 24, dx) for dx in DXS]
                + [(32, 32, dx) for dx in DXS]
                + [(512, 128, 1 / 64)])
    cases = []
    for n in NS:
        ns = sorted({int(round(n * f)) for f in FRACTIONS} | {n})
        for terminal_n in ns:
            for dx in DXS:
                cases.append((n, terminal_n, dx))
    return cases


def atomic_json(path: Path, payload):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def case_path(outdir: Path, spec):
    n, terminal_n, dx = spec
    dx_tag = {1 / 16: '1_16', 1 / 32: '1_32', 1 / 64: '1_64'}[dx]
    return outdir / f'case_N{n}_n{terminal_n}_dx{dx_tag}.json'


def solve_case(spec):
    started = time.monotonic()
    n, terminal_n, dx = int(spec[0]), int(spec[1]), float(spec[2])
    k = n - terminal_n
    if k < 0:
        raise ValueError(f'terminal n exceeds total N: {spec}')
    c0, z, alpha = action_c(n)
    fixed_q = brentq(lambda x: gamma_binomial_cdf_scalar(x, n, c0) - alpha,
                     0.0, MU * n + 12.0 * SIGMA * math.sqrt(n) + 10.0,
                     xtol=2e-11, rtol=2e-14)
    xmax = fixed_q + 8.0 * SIGMA * math.sqrt(n) + 10.0
    count = int(math.ceil(xmax / dx)) + 1
    xmax = dx * (count - 1)
    x = np.arange(count, dtype=np.float64) * dx

    f = gamma_binomial_cdf_vector(x, terminal_n, c0)
    if f[-1] <= alpha:
        raise ArithmeticError(f'terminal CDF does not bracket target for {spec}: F(xmax)={f[-1]}, alpha={alpha}')
    terminal_f = f.copy()
    fixed_dp = terminal_f.copy()
    min_monotonic_increment = math.inf
    min_fixed_dp_increment = math.inf
    max_cdf = float(np.max(f))
    chosen_action_last = None
    controls = np.unique(np.concatenate((ACTION_GRID, np.asarray([c0]))))
    bern = expit(-controls)
    hi_ix, lo_ix, valid_hi, valid_lo, frac = shifted_grid_indices(count, dx, controls)
    one_minus_frac = 1.0 - frac
    c0_index = int(np.argmin(np.abs(controls - c0)))
    if abs(float(controls[c0_index]) - c0) > 2e-14:
        raise AssertionError('fixed action missing from Bellman control list')
    for step in range(k):
        conv = exp_convolution_piecewise_linear(f, dx)
        hi_vals = conv[hi_ix]
        lo_vals = conv[lo_ix]
        shifted = hi_vals * (one_minus_frac[:, None]) + lo_vals * frac[:, None]
        shifted *= valid_hi * one_minus_frac[:, None] + valid_lo * frac[:, None]
        candidates = (1.0 - bern[:, None]) * conv[None, :] + bern[:, None] * shifted
        chosen = np.argmin(candidates, axis=0)
        new_f = candidates[chosen, np.arange(count)]
        new_f[0] = 0.0
        delta_min = float(np.min(np.diff(new_f)))
        min_monotonic_increment = min(min_monotonic_increment, delta_min)
        if delta_min < -2e-12:
            raise ArithmeticError(f'CDF lost monotonicity at step {step}: min increment {delta_min}')
        if not (np.all(np.isfinite(new_f)) and new_f.min() >= -1e-12 and new_f.max() <= 1.0 + 1e-12):
            raise ArithmeticError(f'Bellman CDF left [0,1] at step {step}')
        f = new_f
        fixed_conv = exp_convolution_piecewise_linear(fixed_dp, dx)
        fixed_shift = ((1.0 - frac[c0_index]) * fixed_conv[hi_ix[c0_index]]
                       + frac[c0_index] * fixed_conv[lo_ix[c0_index]])
        fixed_shift *= ((1.0 - frac[c0_index]) * valid_hi[c0_index]
                        + frac[c0_index] * valid_lo[c0_index])
        fixed_dp = ((1.0 - bern[c0_index]) * fixed_conv
                    + bern[c0_index] * fixed_shift)
        fixed_dp[0] = 0.0
        fixed_delta_min = float(np.min(np.diff(fixed_dp)))
        min_fixed_dp_increment = min(min_fixed_dp_increment, fixed_delta_min)
        if fixed_delta_min < -2e-12 or not np.all(np.isfinite(fixed_dp)):
            raise ArithmeticError(f'fixed-action grid reference invalid at step {step}')
        max_cdf = max(max_cdf, float(f[-1]))
        chosen_action_last = chosen

    fixed_grid = gamma_binomial_cdf_vector(x, n, c0)
    q_fixed_grid, fixed_idx, fixed_frac = interp_quantile(x, fixed_grid, alpha)
    q_fixed_dp, fixed_dp_idx, fixed_dp_frac = interp_quantile(x, fixed_dp, alpha)
    q_feedback, feedback_idx, feedback_frac = interp_quantile(x, f, alpha)
    fixed_true_cdf_residual = gamma_binomial_cdf_scalar(fixed_q, n, c0) - alpha
    if k == 0:
        max_grid_difference = float(np.max(np.abs(f - fixed_grid)))
        if max_grid_difference != 0.0:
            raise AssertionError(f'k=0 calibration CDF mismatch: {max_grid_difference}')
    else:
        max_grid_difference = None
    max_feedback_above_fixed_dp = float(np.max(f - fixed_dp))
    if max_feedback_above_fixed_dp > 3e-11:
        raise ArithmeticError(f'feedback CDF exceeds the same-grid fixed-action policy: {max_feedback_above_fixed_dp}')
    if f[-1] <= alpha:
        raise ArithmeticError(f'feedback CDF target not bracketed before xmax: F(xmax)={f[-1]}, alpha={alpha}')

    selected_action = None
    if chosen_action_last is not None:
        action_index = int(chosen_action_last[feedback_idx])
        selected_action = float(controls[action_index])
    row = {
        'status': 'complete',
        'N': n, 'n_terminal': terminal_n, 'k_feedback': k,
        'terminal_fraction': terminal_n / n,
        'grid_spacing': dx, 'grid_points': count, 'x_max': xmax,
        'mu': MU, 'sigma': SIGMA, 'z_N': z, 'target_probability_1_minus_p': alpha,
        'fixed_action_c_N': c0,
        'fixed_action_bernoulli_probability': float(expit(-c0)),
        'action_values_count': int(len(controls)),
        'action_grid_description': '[0,4] step 1/16 plus c_N',
        'terminal_cdf': 'Gamma(n,1)+c_N*Binomial(n, logistic(-c_N))',
        'fixed_reference_cdf': 'Gamma(N,1)+c_N*Binomial(N, logistic(-c_N))',
        'q_fixed_true_gamma_binomial': fixed_q,
        'q_fixed_grid_interpolated': q_fixed_grid,
        'q_fixed_dp_grid_interpolated': q_fixed_dp,
        'q_feedback_grid_interpolated': q_feedback,
        'feedback_gain_vs_true_fixed_quantile': q_feedback - fixed_q,
        'feedback_gain_vs_grid_fixed_quantile': q_feedback - q_fixed_grid,
        'feedback_gain_vs_fixed_dp_quantile': q_feedback - q_fixed_dp,
        'fixed_quantile_grid_interpolation_error': q_fixed_grid - fixed_q,
        'fixed_dp_quantile_error_vs_true_fixed': q_fixed_dp - fixed_q,
        'fixed_true_cdf_quantile_residual': fixed_true_cdf_residual,
        'feedback_cdf_at_xmax': float(f[-1]),
        'target_bracketed': bool(f[0] <= alpha < f[-1]),
        'max_abs_feedback_vs_fixed_grid_cdf': max_grid_difference,
        'max_feedback_cdf_minus_fixed_action_dp_cdf': max_feedback_above_fixed_dp,
        'feedback_cdf_not_above_same_grid_fixed_action_policy': max_feedback_above_fixed_dp <= 3e-11,
        'fixed_grid_quantile_index': fixed_idx,
        'fixed_grid_quantile_fraction': fixed_frac,
        'fixed_dp_quantile_index': fixed_dp_idx,
        'fixed_dp_quantile_fraction': fixed_dp_frac,
        'feedback_grid_quantile_index': feedback_idx,
        'feedback_grid_quantile_fraction': feedback_frac,
        'selected_action_at_feedback_quantile_last_step': selected_action,
        'min_cdf_grid_increment_over_recursion': (min_monotonic_increment if k else 0.0),
        'min_fixed_dp_cdf_grid_increment_over_recursion': (min_fixed_dp_increment if k else 0.0),
        'cdf_min': float(np.min(f)), 'cdf_max': float(np.max(f)),
        'numeric_method': 'float64; analytic Exp convolution of piecewise-linear grid; positive-weight shift interpolation; pointwise minimum over actions at each state',
        'cdf_interpolation_scope': 'state CDF is re-linearized on each uniform grid after each Bellman step; resolution comparison is diagnostic, not a rigorous error enclosure',
        'wall_seconds': time.monotonic() - started,
        'pid': os.getpid(),
    }
    if k == 0:
        row['feedback_gain_vs_grid_fixed_quantile'] = 0.0
        row['feedback_gain_vs_fixed_dp_quantile'] = 0.0
    return row


def solve_and_persist(args):
    spec, outdir_str = args
    row = solve_case(spec)
    outdir = Path(outdir_str)
    outdir.mkdir(parents=True, exist_ok=True)
    atomic_json(case_path(outdir, spec), row)
    return row


def run(outdir: Path, pilot: bool):
    outdir.mkdir(parents=True, exist_ok=True)
    specs = build_cases(pilot)
    start = time.monotonic()
    rows, pids = [], set()
    with cf.ProcessPoolExecutor(max_workers=WORKERS, mp_context=mp.get_context('spawn')) as pool:
        futures = [pool.submit(solve_and_persist, (spec, str(outdir))) for spec in specs]
        for done, future in enumerate(cf.as_completed(futures), start=1):
            row = future.result()
            rows.append(row)
            pids.add(int(row['pid']))
            print(json.dumps({'completed': done, 'total': len(specs), 'N': row['N'],
                              'n_terminal': row['n_terminal'], 'dx': row['grid_spacing'],
                              'gain': row['feedback_gain_vs_true_fixed_quantile'],
                              'seconds': row['wall_seconds']}), flush=True)
    rows.sort(key=lambda r: (r['N'], r['n_terminal'], -r['grid_spacing']))
    if len(rows) != len(specs):
        raise RuntimeError(f'incomplete batch: rows={len(rows)} specs={len(specs)}')
    summaries = {}
    for row in rows:
        key = (row['N'], row['n_terminal'])
        summaries.setdefault(key, {})[str(row['grid_spacing'])] = row
    resolution_pairs = []
    for (n, terminal_n), by_dx in summaries.items():
        ordered_dx = sorted((float(dx) for dx in by_dx), reverse=True)
        for coarse_dx, fine_dx in zip(ordered_dx, ordered_dx[1:]):
            coarse = by_dx[str(coarse_dx)]
            fine = by_dx[str(fine_dx)]
            resolution_pairs.append({
                'N': n, 'n_terminal': terminal_n, 'k_feedback': n-terminal_n,
                'coarse_dx': coarse_dx, 'fine_dx': fine_dx,
                'gain_coarse': coarse['feedback_gain_vs_true_fixed_quantile'],
                'gain_fine': fine['feedback_gain_vs_true_fixed_quantile'],
                'gain_resolution_difference': (fine['feedback_gain_vs_true_fixed_quantile']-
                                               coarse['feedback_gain_vs_true_fixed_quantile']),
            })
    result = {
        'schema': 'finite-prefix-bellman-cdf-comparison-v1',
        'pass_all': all(r['status'] == 'complete' and r['target_bracketed'] for r in rows),
        'mode': 'pilot' if pilot else 'full',
        'scope': 'Finite numerical Gamma-Binomial terminal law and predictable-prefix Bellman comparison only; not a physical strong-maximal norm or theorem.',
        'configured_workers': WORKERS, 'actual_workers': len(pids), 'worker_pids': sorted(pids),
        'case_count': len(rows), 'N_values': sorted(set(r['N'] for r in rows)),
        'terminal_fractions': sorted(set(r['terminal_fraction'] for r in rows)),
        'grid_spacings': list(DXS), 'action_values_count': 66,
        'target_probability': 'alpha_N=sigma*z_N/sqrt(N), where barPhi(z_N)=alpha_N',
        'mu': MU, 'sigma': SIGMA,
        'resolution_pairs': resolution_pairs,
        'cases': rows,
        'wall_seconds': time.monotonic() - start,
        'source_sha256': {str(Path(__file__).resolve()): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        'python_version': platform.python_version(), 'numpy_version': np.__version__,
        'scipy_version': scipy.__version__,
        'resource_note': 'Run under math_guard; RSS and available-memory readings are sampled estimates, not hard OS reservations.',
    }
    summary_path = outdir / ('pilot_results108.json' if pilot else 'results108.json')
    atomic_json(summary_path, result)
    print(json.dumps({'summary': str(summary_path), 'pass_all': result['pass_all'],
                      'cases': len(rows), 'workers': len(pids),
                      'wall_seconds': result['wall_seconds']}, indent=2), flush=True)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--pilot', action='store_true')
    args = ap.parse_args()
    run(args.outdir, args.pilot)


if __name__ == '__main__':
    main()
