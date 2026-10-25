"""Exact symbolic checks for the continuous fixed-radius edge reduction.

The extremal-polytope argument is in proofs/harmonic_cube_minimum.md.
This script checks formulas, not an upper bound on the maximal operator.
"""
import json
import sympy as s

k, t, a, b, c = s.symbols('k t a b c', real=True)
checks = {}
for n in (3, 4):
    y = s.symbols('y0:'+str(n), real=True)
    total = sum(y)
    f = c*total**s.Rational(2-n, 2)-a-b*(sum(v*v for v in y)-3*total**2/s.Integer(n+2))
    directional_hessian = s.diff(f,y[0],2)+s.diff(f,y[1],2)-2*s.diff(f,y[0],y[1])
    assert s.simplify(directional_hessian+4*b) == 0
    edge = c*(k+t)**s.Rational(2-n,2)-a-b*(k+t*t-3*(k+t)**2/s.Integer(n+2))
    expected = b*(6*k-4*t)/5-c/(2*(k+t)**s.Rational(3,2)) if n == 3 else b*(k-t)-c/(k+t)**2
    assert s.simplify(s.diff(edge,t)-expected) == 0
    checks[str(n)] = {'fixed_sum_direction_hessian': '-4*b', 'edge_derivative': str(expected)}
print(json.dumps({'status':'PASS', 'sympy':s.__version__, 'checks':checks,
                  'scope':'symbolic identities only; no new upper bound'},indent=2))
