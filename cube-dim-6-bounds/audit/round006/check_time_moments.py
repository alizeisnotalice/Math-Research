#!/usr/bin/env python3
"""Zero-kernel independent exact/Arb audit of B and time-moment refinement."""
import hashlib
import json
import math
from pathlib import Path
import sympy as sp
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
MAIN = HERE.parents[1]
INPUTS = [MAIN/'coordination/phase_20261025/frozen_thresholds.json',
          MAIN/'upper/round005/cover/010_product_inequality_whole_face.md',
          MAIN/'upper/round006/integral/report.md']
ctx.prec = 256
q = sp.Rational(6, 5)
z = sp.Symbol('z')
F = sum(sp.binomial(63+k, k)*z**k/32**k for k in range(1, 7))
expected = {2:sp.Rational(1167,200),4:sp.Rational(187576101,2560000),6:sp.Rational(73354040202609,16384000000)}
eps = arb(1)/1000
B = arb(1)/16000
A = arb(243)**(-arb(1)/6)*arb.pi()**(-arb(11)/6)
K = arb(6)/7*A
rows=[]
for m in (2,4,6):
    poly=sum(abs(sp.binomial(q,j))*F**j for j in range(1,m+1))
    D=sp.factorial(m)*sp.Poly(poly,z).coeff_monomial(z**m)
    assert D == expected[m]
    rho=D/(2**m*sp.factorial(m))
    dar=arb(int(D.p))/int(D.q)
    rar=arb(int(rho.p))/int(rho.q)
    a=arb(m)-arb(7)/6
    # J_m assigned envelope: integral_eps^infinity A*t^(1/6-m)/(1+1024*t^2) dt.
    # Exact analytic lower from [eps,2eps], upper by dropping denominator.
    integral_upper=A*eps**(-a)/a
    integral_lower=integral_upper*(1-arb(2)**(-a))/(1+4096*eps**2)
    ratio_lower=integral_lower*eps**m/B
    ratio_upper=integral_upper*eps**m/B
    # Leading one-axis Peano error / target I=B^q with h=2eps.
    theta_low=rar*arb(2)**m*ratio_lower
    theta_up=rar*arb(2)**m*ratio_upper
    # Bounds on h/eps for leading allocated error <=10% target.
    h_ratio_lower=(arb(1)/10/(rar*ratio_upper))**(arb(1)/m)
    h_ratio_upper=(arb(1)/10/(rar*ratio_lower))**(arb(1)/m)
    rows.append(dict(m=m,D_exact=str(D),rho_exact=str(rho),
      scalar_moment_budget_lower=integral_lower.str(45),scalar_moment_budget_upper=integral_upper.str(45),
      improvement_factor_lower=ratio_lower.str(45),improvement_factor_upper=ratio_upper.str(45),
      theta_at_h_2eps_lower=theta_low.str(45),theta_at_h_2eps_upper=theta_up.str(45),
      h_over_eps_for_10pct_leading_budget_lower=h_ratio_lower.str(45),
      h_over_eps_for_10pct_leading_budget_upper=h_ratio_upper.str(45)))

# Correct algebraic sharp pointwise product bound for coordinate second derivative.
r=sp.Symbol('r', nonnegative=True)
phi=(65*r-1)/(16*(1+r)**2)
critical=sp.Rational(67,65)
assert sp.diff(phi,r).subs(r,critical)==0
assert sp.factor(phi.subs(r,critical))==sp.Rational(4225,4224)
output=dict(status='PASS_EXACT_HIGH_DERIVATIVE_AND_TIME_MOMENT_BUDGET_AUDIT',
 scope='No kernel evaluations. Scalar envelope only; no certified I, g, or C6 improvement.',
 epsilon='1/1000',norm_reference='1/16000 (conditional target, not an achieved norm)',
 tail_budget=(K*eps**(arb(7)/6)).str(45),rows=rows,
 sharp_product_second_derivative='4225/4224',
 inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUTS},
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'check_time_moments_results.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
