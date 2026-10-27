#!/usr/bin/env python3
"""Exact rational spot checks for the independent R030 proof.

The general theorem is proved analytically in proof.md. These finite checks
verify algebra and threshold arithmetic; they do not prove a continuum claim.
"""
from fractions import Fraction as Q
from math import comb, factorial
import hashlib
import json
from pathlib import Path
import platform
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def integral_weight(n, a, moment):
    return sum((Q((-1)**k * comb(n-1, k), n**k)
                * a**(k+moment+1) / (k+moment+1)
                for k in range(n)), Q(0))

assert sum((Q(3**k, factorial(k)) for k in range(6)), Q(0)) == Q(92, 5)
assert Q(92, 5) > 18
assert 16 * Q(3, 4)**16 < 1
assert Q(3, 4) * Q(17, 16) < 1
assert 18**3 - 120*18 - 240 == 3432
assert 3*18**2 - 120 > 0
assert Q(1, 2) - Q(3, 18) > 0

cases = []
for n in (1, 2, 7, 16):
    assert integral_weight(n, Q(n), 0) == 1
    assert integral_weight(n, Q(n), 1) == Q(n, n+1)
    for x in (Q(1, 3), Q(n, 2), Q(3*n, 4)):
        a = (n-x)/2
        h = (1+x/n)/2
        assert integral_weight(n, a, 0) == 1-h**n
        assert 0 < integral_weight(n, a, 1) < 1
        for u in (Q(0), a/3, a):
            rho = Q(n)/(x+2*u)
            ell = 1+Q(x, 2*n)-Q(1, 2*rho)
            assert ell == 1-u/n
            # -rho'(u)/rho(u) = 1/(u+x/2).
            assert Q(2, x+2*u) == Q(1)/(u+x/2)
        cases.append({'n': n, 'x': str(x)})

# The derivative of y+y^2+log(1-y) is y(1-2y)/(1-y).
for y in (Q(0), Q(1, 7), Q(1, 3), Q(1, 2)):
    assert 1+2*y-Q(1)/(1-y) == y*(1-2*y)/(1-y)
    assert y*(1-2*y)/(1-y) >= 0

source_paths = [
    'AGENTS.md',
    'cube-dim-7-bounds/revision_v3/cover/audit.md',
    'cube-dim-7-bounds/revision_v3/cover/asymptotic_cross_review.md',
    'cube-dim-7-bounds/revision_v3/asymptotic/proof.md',
    'cube-dim-7-bounds/audit/round030/uniform/proposal_snapshot.md',
]
sources = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
           for p in source_paths}
result = {
    'status': 'PASS', 'python': sys.version, 'platform': platform.platform(),
    'arithmetic': 'Python stdlib fractions.Fraction (exact)',
    'base_exp3_partial_sum': '92/5',
    'n16_tail_comparison': str(16*Q(3,4)**16),
    'absorption_polynomial_at_18': 3432,
    'finite_integral_cases': cases,
    'scope': 'finite exact checks only; general statements rely on proof.md',
    'source_sha256': sources,
}
(HERE/'check_result.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
