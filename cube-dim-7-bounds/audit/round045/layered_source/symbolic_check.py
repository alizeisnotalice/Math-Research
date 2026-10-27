#!/usr/bin/env python3
"""New symbolic identities and exact constants, no parameter integrations."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform
import sys

sys.path.insert(0, '/Users/zhengzhihao/Documents/ChatGPT/bnupt/tmp/y_pendulum_full/deps')
import sympy as s

HERE = Path(__file__).resolve().parent
SOURCE = Path('/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round045/layered_source/PROOF.md')
EXPECTED = '22db93b275416b90e9f15a0118a88cd412fee0c2d0b20c4cce1ed32769f259ac'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED

# The tail variable y is |Z|-(rho-g), with 0<=y<=g< R.
R, rho, g, y, b = s.symbols('R rho g y b', positive=True)
loss_mean = s.integrate(y/(2*R), (y, 0, g))/rho
assert s.simplify(loss_mean - g**2/(4*R*rho)) == 0

# Inverse-coordinate mean F and excluded-tail mass E(b), 0<=b<=g.
# Verify by derivatives and boundary conditions, avoiding log branch rewrites.
Fminus1 = (2*R*s.log(2*R/(2*R-g))-g)/rho
tail_integrand = (1/(1-y/(2*R))-1)/rho
assert s.simplify(s.diff(Fminus1,g)-tail_integrand.subs(y,g)) == 0
assert s.simplify(Fminus1.subs(g,0)) == 0
E = R/rho*s.log((2*R-g+b)/(2*R-g))
assert s.simplify(s.diff(E,b) - 1/(2*rho*(1-(g-b)/(2*R)))) == 0
assert s.simplify(E.subs(b,0)) == 0

# Primitive of E(b): the average amount removed from the reference J row
# is exactly F-1, so the reference mean is exactly one.
v = 2*R-g+b
primitive = R/rho*(v*s.log(v/(2*R-g))-v)
assert s.simplify(s.diff(primitive,b)-E) == 0
reference_removed = (primitive.subs(b,g)-primitive.subs(b,0))/R
assert s.simplify(reference_removed-Fminus1) == 0

# Full finite-cover objective at its chosen tau, checked algebraically.
A, C, q = s.symbols('A C q', positive=True)
full_cost = (1+A/q)*(q+C)
assert s.simplify(full_cost.subs(q,s.sqrt(A*C))-(s.sqrt(A)+s.sqrt(C))**2) == 0

d, eps, K = F(1,2**64), F(1,2**10), F(2**20)
B, HR = 4*K, F(2**24)
remaining = F(1,64)-3*F(1,4096)-F(1,131072)
assert remaining >= F(1,128)
assert 8/K == F(1,2**17)
assert F(4096)/HR == F(1,2**12)
assert (F(1,128)-F(1,256))/HR == F(1,2**32)
assert 2*B*d*d == (F(1,2**32)*d)/512
assert F(1,2**32)*d/2 == F(1,2**97)
assert B/2-K == K > 0
assert B*d*d < F(1,2)
result = {
    'proof_sha256': EXPECTED,
    'python': platform.python_version(), 'sympy': s.__version__,
    'method': 'SymPy symbolic derivative/boundary identities and Fraction constants; no integration samples or parameter scan',
    'identities': ['free-coordinate mean loss', 'inverse-coordinate mean F', 'one-sided excluded mass', 'reference E[J]=1', 'complete cover cost at optimal tau'],
    'remaining_energy_coefficient': str(remaining),
    'retained_energy_lower_bound': str(F(1,128)),
    'deficit_coefficient_over_nd': str(F(1,2**32)),
    'scalar_buffer_over_deficit': '1/512',
    'final_gain_over_n': str(F(1,2**97)),
    'all_assertions_passed': True,
}
(HERE/'symbolic_result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
