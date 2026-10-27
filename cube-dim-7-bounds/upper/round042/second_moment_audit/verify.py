"""Exact one-coordinate and rational checks; no production numerics."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sympy as s

root=Path(__file__).resolve().parent
v,D,r,rho=s.symbols('v D r rho',positive=True)
ex=s.integrate((D-v)/rho,(v,0,D))
ex2=s.integrate((D-v)**2/rho,(v,0,D))
assert s.factor(ex-D**2/(2*rho))==0
assert s.factor(ex2-D**3/(3*rho))==0
ell2=s.integrate((1-(D-v)/(2*r))**2/rho,(v,0,D))+(rho-D)/rho
assert s.factor(ell2-(1-D**2/(2*r*rho)+D**3/(12*r*r*rho)))==0
checks=['exact free-coordinate first and second deficit moments',
        'exact E ell^2 with positive cubic term']

d=s.symbols('d',positive=True)
q=1-d/4
D0=d-d*d/100
r0=q+D0
x=d**-3-d**-1
edge=s.limit(-x*D0**2/(2*r0*q)+x*D0**3/(12*r0*r0*q)+1/(2*d),d,0,dir='+')
assert edge==s.Rational(103,300)
checks.append('edge constant 103/300 at q=q_max')

assert 6400<15*4**13
assert F(1,8)+F(13,50)+F(1,1024)<F(2,5)
assert F(2,5)+F(4,39)<1
assert F(8,13)<1
assert F(19,16)+F(1,800)<2
assert 1-F(1,4)/2==F(7,8)>F(1,2)
assert 16*3==48
checks.append('all uniform L>=64 rational gates and C=48')
receipt={
 'status':'PASS: exact moment and uniform-envelope algebra',
 'condition':'log(n)>=64',
 'second_moment_constant':48,
 'checks':checks,
 'sympy':s.__version__,
 'scope':'moment/column-upper audit only; no new operator bound claimed',
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                  for name in ('PROOF.md','verify.py')},
 'no_production_search':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
