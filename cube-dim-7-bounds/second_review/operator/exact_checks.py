#!/usr/bin/env python3
"""Bounded exact algebra/record checks only; no transcendental evaluation or search."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json, sys, time

root = Path(__file__).resolve().parents[2]
start = time.perf_counter()
expected = [Q(-3), Q(15,8), Q(-5,6), Q(15,64), Q(-3,80), Q(1,384)]
actual = [Q(comb(6,j), j) * Q(-1,2)**j for j in range(1,7)]
assert actual == expected
r = json.loads((root / 'upper/round028/cap_candidate/result.json').read_text())
t = Q(3391749,524288)
q, h = 1 + 1/(2*t), (1 + 1/t)/2
coef = 7-Q(7,2)*q**6
const = 1-h**7-Q(7,2)*sum(
    Q(comb(6,j),j)*q**(6-j)*Q(-1,2)**j*(1-t**(-j)) for j in range(1,7))
assert Q(r['exact_t']) == t
assert Q(r['log_coefficient']) == coef
assert Q(r['rational_constant']) == const
ep = r['candidate_interval']['upper']
upper = Q(int(ep[0]))*Q(2)**int(ep[1])
decimal = Q(24554960525965,10**12)
assert upper < decimal
assert Q(r['candidate_grid_upper']) == decimal
out = dict(
    status='PASS', python=sys.version, wall_seconds=time.perf_counter()-start,
    exact_binomial_coefficient_checks=6, exact_registered_t_rows=1,
    checks=['six explicit antiderivative coefficients',
            'recorded R028 rational log coefficient and constant',
            'recorded dyadic upper endpoint strictly below stated decimal'],
    log_evaluations=0, sqrt_evaluations=0, integral_evaluations=0,
    optimization_runs=0, parameter_scans=0, new_candidate_points=0,
    limitation='Does not independently enclose log/sqrt or certify the stored interval.')
print(json.dumps(out,indent=2))
