#!/usr/bin/env python3
"""One preregistered rational SPD matrix, one fixed sign flip, one positive series."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import argparse
import hashlib
import json
import os
import time
import sympy as sp
import flint
from flint import arb, ctx

BASE = Path(__file__).resolve().parent
START = time.monotonic()
parser = argparse.ArgumentParser()
parser.add_argument('--output-dir', type=Path, default=BASE)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(name, value):
    (args.output_dir/name).write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def checkpoint(stage):
    elapsed = time.monotonic()-START
    save('checkpoint.json', {'pid': os.getpid(), 'stage': stage,
                            'elapsed_seconds': elapsed, 'limit_seconds': 120,
                            'resume': 'python -B verify.py'})
    if elapsed > 120:
        raise TimeoutError('Recoverable checkpoint retained')


def ar(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def ends(x):
    return F(str(x.lower().fmpq())), F(str(x.upper().fmpq()))


def decimal(x, up=False, digits=55):
    scale = 10**digits
    value = x*scale
    n = value.__ceil__() if up else value.__floor__()
    sign = '-' if n < 0 else ''
    integer, fraction = divmod(abs(n), scale)
    return f'{sign}{integer}.{fraction:0{digits}d}'


def interval(lo, hi):
    assert lo <= hi
    return {'lower_fraction': str(lo), 'upper_fraction': str(hi),
            'decimal_outward': [decimal(lo), decimal(hi, True)]}


def root_bounds(x, degree, digits=110):
    x = F(x)
    scale = 10**digits
    target = x.numerator*scale**degree//x.denominator
    lo, hi = 0, 1 << ((target.bit_length()+degree-1)//degree)
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**degree <= target:
            lo = mid
        else:
            hi = mid
    lower = F(lo, scale)
    upper = lower if lower**degree == x else F(lo+1, scale)
    assert lower**degree <= x <= upper**degree
    return lower, upper


checkpoint('registered_start')
registration = json.loads((BASE/'registration.json').read_text())
source = registration['sources'][0]
assert sha(Path(source['path'])) == source['sha256']
q = F(5, 4)
M = sp.Rational(3, 4)*sp.eye(5)+sp.Rational(1, 4)*sp.ones(5)
S = sp.diag(-1, 1, 1, 1, 1)
A = M.inv()
assert A == sp.Rational(4, 3)*sp.eye(5)-sp.Rational(1, 6)*sp.ones(5)
assert M.det() == sp.Rational(81, 128)
assert M.eigenvals() == {sp.Rational(3, 4): 4, sp.Integer(2): 1}
assert all(M[i, i] == 1 for i in range(5))
assert S*M*S != M
body_vectors = [sp.Matrix((1,)+e) for e in product((-1, 1), repeat=4)]
body_weights = [sp.Rational(19, 64) if tuple(v) == (1,)*5 else sp.Rational(3, 64)
                for v in body_vectors]
assert sum(body_weights) == 1
reconstructed = sp.zeros(5)
for w, v in zip(body_weights, body_vectors):
    reconstructed += w*v*v.T
assert reconstructed == M
t = sp.symbols('t', real=True)
det_ratio = sp.factor(M.det()*(t*A+(1-t)*S*A*S).det())
assert sp.simplify(det_ratio-(1+sp.Rational(2, 3)*t*(1-t))) == 0
kappa = F(str(A[0, 0]-1))
s = kappa/(1+kappa)
assert (kappa, s) == (F(1, 6), F(1, 7))
assert F(5, 14)*s == F(5, 98)
checkpoint('exact_matrix_identities_complete')

# One positive series, fixed truncation. No quadrature or optimizer.
N = 48
term, partial = F(1), F(1)
for n in range(N):
    ratio = s*(F(n)+F(1, 2))**2/((F(n)+F(7, 4))*(n+1))
    assert 0 < ratio < s
    term *= ratio
    partial += term
next_ratio = s*(F(N)+F(1, 2))**2/((F(N)+F(7, 4))*(N+1))
next_term = term*next_ratio
tail = next_term/(1-s)
assert 0 < tail < F(1, 10**43)
ctx.prec = 384
sqrt_ball = ar(1-s).sqrt()
rho_lo = ends(sqrt_ball*ar(partial))[0]
rho_hi = ends(sqrt_ball*ar(partial+tail))[1]
assert 0 < rho_lo < rho_hi < 1
gap_lo, gap_hi = 1-rho_hi, 1-rho_lo
assert gap_lo > F(5, 98)
beta_lo, beta_hi = ends(1-ar(F(1, 2)).sqrt().sqrt())
assert 0 < beta_lo <= beta_hi < 1
loss_lo, loss_hi = beta_lo*gap_lo, beta_hi*gap_hi
factor_lo, factor_hi = 1-loss_hi, 1-loss_lo
norm_factor_lo = root_bounds(factor_lo**4, 5)[0]
norm_factor_hi = root_bounds(factor_hi**4, 5)[1]
assert 0 < norm_factor_lo <= norm_factor_hi < 1
pair_energy_lo, pair_energy_hi = (1+rho_lo)/2, (1+rho_hi)/2
holder_factor_lo = root_bounds(pair_energy_lo, 4)[0]
holder_factor_hi = root_bounds(pair_energy_hi, 4)[1]
assert factor_hi < holder_factor_lo  # sharper pair inequality versus generic Holder.
checkpoint('one_dimensional_positive_series_complete')

# The sphere normalization is local to this same matrix, not a new model test.
d = F(81, 128)
omega = ar(F(8, 3))*arb.pi()**2
norm_q_ball = omega*ar(1/d).sqrt().sqrt().sqrt()
A_lo, A_hi = ends(norm_q_ball)
B_lo, B_hi = A_lo*rho_lo, A_hi*rho_hi
face_prefactor_lo, face_prefactor_hi = A_lo/160, A_hi/160
results = {
    'status': 'single_preregistered_matrix_and_flip_verified',
    'matrix_count': 1, 'sign_flip_count': 1,
    'dimension': 5, 'q': '5/4', 'det_M': str(d),
    'M_eigenvalues': {'3/4': 4, '2': 1},
    'M_inverse': '(4/3)I-(1/6)11^T',
    'body_shape_admissibility': {'all_plus_weight': '19/64', 'other_fifteen_weights': '3/64'},
    'determinant_ratio': str(det_ratio), 'kappa': str(kappa), 's': str(s),
    'positive_series': {'first_index': 0, 'last_index': N,
                        'partial_sum_fraction': str(partial),
                        'first_omitted_term_fraction': str(next_term),
                        'tail_upper_fraction': str(tail),
                        'all_successive_ratios_below': str(s)},
    'rho_normalized_cross_energy': interval(rho_lo, rho_hi),
    'energy_gap_one_minus_rho': interval(gap_lo, gap_hi),
    'universal_rational_gap_lower': '5/98',
    'beta_q_one_minus_two_power_one_minus_q': interval(beta_lo, beta_hi),
    'pair_q_power_relative_loss_expression': interval(loss_lo, loss_hi),
    'pair_q_power_relative_upper_factor_expression': interval(factor_lo, factor_hi),
    'pair_norm_relative_upper_factor_expression': interval(norm_factor_lo, norm_factor_hi),
    'two_member_holder_q_power_factor_expression': interval(holder_factor_lo, holder_factor_hi),
    'sphere_integral_H_to_q': interval(A_lo, A_hi),
    'sphere_cross_integral': interval(B_lo, B_hi),
    'face_full_orbit_unimproved_q_prefactor': interval(face_prefactor_lo, face_prefactor_hi),
    'global_R_expectation_computed': False,
    'full_group_integral_computed': False,
    'new_numeric_I_bound': False, 'new_C5_bound': False,
    'note': 'Norm factors are rigorously enclosed upper-bound expressions, not measured actual pair norms.'
}
save('results.json', results)
checkpoint('complete')
receipt = {'pid': os.getpid(), 'elapsed_seconds': time.monotonic()-START,
           'exit_status': 0, 'precision_bits': ctx.prec,
           'python_flint': flint.__version__, 'sympy': sp.__version__,
           'results_sha256': sha(args.output_dir/'results.json'),
           'script_sha256': sha(Path(__file__)), 'registration_sha256': sha(BASE/'registration.json')}
save('run_receipt.json', receipt)
print(json.dumps(receipt, sort_keys=True))
print(json.dumps({'rho': results['rho_normalized_cross_energy']['decimal_outward'],
                  'pair_q_power_loss': results['pair_q_power_relative_loss_expression']['decimal_outward'],
                  'pair_norm_upper_factor': results['pair_norm_relative_upper_factor_expression']['decimal_outward']},
                 sort_keys=True))
