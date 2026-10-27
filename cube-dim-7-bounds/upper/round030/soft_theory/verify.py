#!/usr/bin/env python3
"""Small algebra/quadrature checks; no parameter optimization."""
import hashlib
import json
import math
import platform
from pathlib import Path

import mpmath as mp
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
mp.mp.dps = 40
h, beta, v, z = s.symbols('h beta v z', positive=True)
checks = []
for n in range(1, 9):
    power = sum(s.binomial(n, k)*beta*h**k/(beta+k) for k in range(n+1))
    independent_power = sum(s.binomial(n, k)*h**k*s.integrate(beta*v**(beta+k-1), (v, 0, 1)) for k in range(n+1))
    exponential = sum(s.binomial(n, k)*s.factorial(k)*h**k for k in range(n+1))
    independent_exponential = sum(s.binomial(n, k)*h**k*s.gamma(k+1) for k in range(n+1))
    checks.append({'n': n, 'power_coverage': s.simplify(power-independent_power)==0,
                   'exponential_coverage': s.simplify(exponential-independent_exponential)==0})

# The difference between AB and (E sqrt(a b))^2 for a two-point law.
p, x1, x2, y1, y2 = s.symbols('p x1 x2 y1 y2', real=True)
gap = (p*x1*x1+(1-p)*x2*x2)*(p*y1*y1+(1-p)*y2*y2)-(p*x1*y1+(1-p)*x2*y2)**2
cs_exact = s.factor(gap-p*(1-p)*(x1*y2-x2*y1)**2)==0

def b_from_v(n, vv):
    """Hard B_n(t), v=1/t in (0,1); q>=2 has residual zero."""
    vv = mp.mpf(vv)
    if vv >= 1:
        return mp.mpf(0)
    if vv <= 0:
        return mp.inf
    c = 1+vv/2
    integral = c**(n-1)*(-mp.log(vv))
    for k in range(1, n):
        integral += math.comb(n-1,k)*c**(n-1-k)*(-mp.mpf(1)/2)**k*(1-vv**k)/k
    return 1-((1+vv)/2)**n-n*mp.log(vv)-n*integral/2

def budget(a, b):
    return (mp.sqrt(a)+mp.sqrt(b))**2

def out(x):
    return mp.nstr(x, 32)

n = 7
q_values = [mp.mpf('1.12'), mp.mpf('1.20')]
aa = [q**n for q in q_values]
bb = [b_from_v(n,q-1) for q in q_values]
ff = [budget(a,b) for a,b in zip(aa,bb)]
mixed = budget(sum(aa)/2, sum(bb)/2)
two_point = {'q': list(map(out,q_values)), 'component_F': list(map(out,ff)),
             'mixed_F': out(mixed), 'average_component_F': out(sum(ff)/2),
             'CS_gap': out(mixed-sum(ff)/2)}

hh, be = mp.mpf('0.2'), mp.mpf('2')
ap = sum(math.comb(n,k)*be*hh**k/(be+k) for k in range(n+1))
bp = mp.quad(lambda w: be*w**(be-1)*b_from_v(n,hh*w), [0,mp.mpf('0.1'),1])
ae = sum(math.comb(n,k)*math.factorial(k)*hh**k for k in range(n+1))
be_cost = mp.quad(lambda w: mp.exp(-w/hh)/hh*b_from_v(n,w), [0,mp.mpf('0.1'),1])
continuous = {'n':n, 'h':out(hh), 'power_beta':out(be),
              'power_A':out(ap), 'power_B':out(bp), 'power_F':out(budget(ap,bp)),
              'exponential_A':out(ae), 'exponential_B':out(be_cost),
              'exponential_F':out(budget(ae,be_cost))}

# Independently integrate the original rho formula at the two hard samples.
hard_integral_checks = []
for q in q_values:
    t = 1/(q-1)
    ii = mp.quad(lambda rho:(1+1/(2*t)-1/(2*rho))**(n-1)/rho,[1,t])
    original = 1-((1+1/t)/2)**n+n*mp.log(t)-n*ii/2
    err = abs(original-b_from_v(n,q-1))
    hard_integral_checks.append({'q':out(q),'absolute_error':out(err),'passed':err<mp.mpf('1e-35')})

sources = [HERE/'proof.md', HERE/'verify.py', ROOT/'paper/main.tex', ROOT/'revision_v3/asymptotic/proof.md']
receipt = {'status':'passed', 'python':platform.python_version(), 'sympy':s.__version__,
           'mpmath':mp.__version__, 'mpmath_dps':mp.mp.dps,
           'source_checkout':'f65a5bf92911ebed4cf803071bd7693978e5e582',
           'coverage_checks':checks, 'two_point_CS_identity':cs_exact,
           'hard_integral_checks':hard_integral_checks,
           'two_point_example':two_point, 'continuous_examples':continuous,
           'noncommon_corner_counterexample':{
               'n':1, 'a':10, 'q':'11/10', 'alpha_and_c':'10/11', 't':10,
               'H_origin':out(mp.mpf(1)/11+mp.log(10)),
               'H_corner':out(mp.mpf(1)/2+mp.log(10)/2),
               'origin_minus_corner':out(mp.log(10)/2-mp.mpf(9)/22),
               'exact_positive_lower_bound':'1/11, using log(10)>1'},
           'interpretation':'Small symbolic and numerical checks corroborate formulas; they do not replace the analytic all-n proof and perform no optimization.',
           'sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}}
assert cs_exact and all(all(row[k] for k in ('power_coverage','exponential_coverage')) for row in checks)
assert all(row['passed'] for row in hard_integral_checks)
assert mixed >= sum(ff)/2
(HERE/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'two_point':two_point,'continuous':continuous},indent=2))
