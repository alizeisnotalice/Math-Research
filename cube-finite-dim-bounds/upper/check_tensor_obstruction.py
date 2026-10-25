"""Finite U01 certificate: rational negative hyperplane masses + Fourier check."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sympy as sp

x = sp.symbols('x', real=True)
integral = sp.integrate(1-x/2, (x, sp.Rational(1,2), 1))
assert integral == sp.Rational(5,16)
rows = []
for n in (3,4):
    negative_mass = -Fraction(5,16)**(n-1)
    # E=(-1/4,1/4) x (1/2,1)^(n-1), excludes zero and ±2 planes.
    # k''=(δ_{-2}+δ_2)/2-δ_0, hence only j=1 contributes.
    assert negative_mass < 0
    rows.append({'n':n, 'tent_mass':2**n,
                 'laplacian_slab_mass_exact':str(negative_mass),
                 'laplacian_slab_mass_decimal':float(negative_mass)})

# Normalize each tent to probability. q(t)=(sin(t)/t)^2.
# F(xi)=sum(xi_j^2) product q(xi_j); any coordinate-reflection
# invariant Levy exponent obeys |F(s,t)-F(s,0)| <= F(0,t).
s,t=sp.pi/2,sp.pi
q=lambda z: sp.simplify((sp.sin(z)/z)**2) if z!=0 else sp.Integer(1)
f_s0=sp.simplify(s*s*q(s))
f_0t=sp.simplify(t*t*q(t))
f_st=sp.simplify((s*s+t*t)*q(s)*q(t))
assert (f_s0,f_0t,f_st)==(1,0,0)
distance=sp.Matrix([[0,f_s0,f_0t,f_st],
                    [f_s0,0,f_st,f_0t],
                    [f_0t,f_st,0,f_s0],
                    [f_st,f_0t,f_s0,0]])
coeff=sp.Matrix([1,1,-1,-1])
cnd_defect=(coeff.T*distance*coeff)[0]
assert sum(coeff)==0 and cnd_defect==4
result={'status':'strict counterexample certificate', 'rows':rows,
        'fourier_exact':{'s':'pi/2','t':'pi','F_s0':str(f_s0),
                         'F_0t':str(f_0t),'F_st':str(f_st),
                         'required_inequality':'1 <= 0 (false)',
                         'four_point_coefficients':list(coeff),
                         'four_point_cnd_quadratic_form':str(cnd_defect)},
        'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(result,indent=2,default=str))
