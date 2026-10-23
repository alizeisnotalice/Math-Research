"""Algebraic specialization of the source identities; assumes R>1, 0<r<1."""
import sympy as s
import json
from pathlib import Path
r,R=s.symbols('r R', positive=True)
L=s.symbols('log_one_minus_r_over_R')
tail=2*(2*(-R*L/r**2-1/r)+1/R)/(R-1)
j=s.Rational(1,2)-R/r-R**2*L/r**2
assert s.factor(tail-4*j/(R*(R-1)))==0
# Log identities valid on 0<r<1; use independent symbols after expansion.
lr,lc,l2,ld=s.symbols('log_r log_one_minus_r log_two log_two_minus_r')
t=r/(2-r)
f=2-(1-t)*(l2+lc-ld)-(1+t)*(l2-ld)+2*(lr-ld)
T=2*(2-r)*(1-l2+lr)-2*(1-r)*lc
assert s.factor((2-r)*f-T)==0
out={'tail_identity_residual':'0','bump_identity_residual':'0',
     'assumptions':'R>1, 0<r<1; real logarithm product identities',
     'status':'symbolic identities verified; continuity argument is in report'}
Path(__file__).with_name('symbolic_diagonal_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
