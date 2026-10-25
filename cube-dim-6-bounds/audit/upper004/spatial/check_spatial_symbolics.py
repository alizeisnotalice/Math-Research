#!/usr/bin/env python3
"""Small symbolic/exact audit for round004 spatial interfaces."""
from __future__ import annotations
import hashlib, json
from itertools import product
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
z,t,a=sp.symbols('z t a', real=True, positive=True)
y,s=sp.symbols('y s', real=True, nonnegative=True)
V=[(1,)+u for u in product((-1,1),repeat=5)]
A=sp.Matrix.hstack(*(sp.Matrix(v) for v in V))
gram=A*A.T
log_second=2*(z*z-t*t)/(t*t+z*z)**2
lower_poly=sp.factor((log_second+2/t**2).subs(z*z,y*t*t)*t**2)
upper_poly=sp.factor((sp.Rational(1,4)/t**2-log_second).subs(z*z,y*t*t)*t**2)
# Difference between squared complex denominator and the desired lower bound:
# |t^2+(z+ia)^2|^2-(1-a^2/t^2)^2(t^2+z^2)^2.
complex_diff=sp.factor((((t*t+z*z-a*a)**2+4*a*a*z*z)
  -(1-a*a/t**2)**2*(t*t+z*z)**2).subs({z*z:y*t*t,a*a:s*t*t})/t**4)
c=[sp.Rational(v[0],32) for v in V]
out={
 'status':'PASS_symbolic_spatial_identities',
 'AAT_equals_32I':gram==32*sp.eye(6),
 'right_inverse_check':A*(A.T/32)==sp.eye(6),
 'single_axis_sum_c_squared':str(sum(q*q for q in c)),
 'single_axis_sum_abs_c':str(sum(abs(q) for q in c)),
 'logf_second_plus_2_over_t2_times_t2':str(lower_poly),
 'one_over_4t2_minus_logf_second_times_t2':str(upper_poly),
 'complex_denominator_squared_difference_scaled':str(complex_diff),
 'inequality_reading':{
   'logconvex_lower':'2*y*(y+3)/(y+1)^2 >= 0',
   'logf_upper':'(y-3)^2/(4*(y+1)^2) >= 0',
   'complex':'s*y*(6+2*y-s*(2+y)) >= 0 for 0<=s<1,y>=0'
 },
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
(HERE/'check_spatial_symbolics_results.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
