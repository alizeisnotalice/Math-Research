#!/usr/bin/env python3
"""Small exact constants only; no dimension grid, optimization or MATLAB."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib
import json
import platform
import sys

y = 64
assert Q(81, 127) < Q(2, 3)
assert 8*Q(9, 16*y) <= Q(1, 8)
assert Q(5, 16)-Q(1, 100) > Q(1, 4)
exp_partial = sum((Q(9, 2)**k/factorial(k) for k in range(6)), Q(0))
assert exp_partial > 50
assert 4*y**3 >= 192
assert Q(1, 4*y*(1-Q(1, 2*y))) < 1
large_base = 256*(Q(91, 32768)+Q(8101, 32768**3))
assert large_base < 1
# Only moderate rational parameters test the shell-face identities.
for n in (4, 16, 64):
    t = Q(n, 2)
    beta, zeta = 1+1/t, Q(1, n**3)
    kappa = beta-zeta
    assert zeta <= (beta-1)/2 and kappa > 1
    a = 8*n**4
    alpha = Q(a, a+1)
    assert alpha*beta-kappa-Q(2, a+1) > 0
    for rho in (Q(4), Q(7)):
        finite = alpha-(alpha*beta-Q(1, a+1))*rho
        limiting = 1-beta*rho
        assert finite-limiting == ((beta+1)*rho-1)/(a+1) > 0
    for w in (Q(-1), Q(0), Q(1)):
        rho = Q(4)
        side_beta = min(Q(1), (1+beta)/2-abs(w)/(2*rho))
        side_kappa = min(Q(1), (1+kappa)/2-abs(w)/(2*rho))
        assert 0 <= side_beta-side_kappa <= zeta/2

here = Path(__file__).resolve().parent
result = {
    'status': 'PASS', 'arithmetic': 'stdlib fractions.Fraction, exact',
    'python': sys.version, 'platform': platform.platform(),
    'exp_9_over_2_partial_sum_k0_to5': str(exp_partial),
    'large_L_base_rational_upper_bound': str(large_base),
    'scope': 'small exact constants and rational identities only; no search',
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(here/'constant_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
