#!/usr/bin/env python3
"""Small exact checks for R002; the general results are proved in research.md."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json
import sympy as s

x,L,D=s.symbols('x L D', positive=True)
polynomials=[]
for n in range(1,9):
    # R(x)-R(0), after replacing D by its upper bound n L/2.
    lower=s.expand((1+n*L)*((1+x)**n-1)-s.Rational(n,2)*L*x+x**n)
    p=s.Poly(lower,x,L)
    assert all(v>0 for v in p.coeffs())
    assert lower.subs(x,0)==0
    polynomials.append({'n':n,'lower_difference':str(lower),
                        'all_nonzero_coefficients_positive':True})
a,t,r=s.symbols('a t r', positive=True)
alpha=a/(a+1); beta=1/(a+1); q=1+1/t; c=alpha*q-beta
identities={
 'center_threshold':s.factor(alpha-alpha*q+alpha/t),
 'cap_threshold':s.factor(alpha*q-beta-c),
 'c_minus_one':s.factor(c-1-(a/t-2)/(a+1))}
assert all(v==0 for v in identities.values())

a0=Q(20);t0=Q(3);al=a0/(a0+1);be=1/(a0+1);qq=1+1/t0;cc=al*qq-be
checked=0
for n in (1,2):
 for rr in (Q(1),Q(3,2),Q(5,2)):
  us=[al*Q(j,4) for j in range(5)]
  zs=[rr*j for j in (Q(-1),Q(-3,4),Q(-1,2),Q(0),Q(1,2),Q(3,4),Q(1))]
  for u in product(us,repeat=n):
   for z in product(zs,repeat=n):
    lo=[max(-al*rr,zi-be*rr) for zi in z]
    hi=[min(al*rr,zi+be*rr) for zi in z]
    intervals=all(l<h for l,h in zip(lo,hi))
    original=intervals and any(l<ui-al*qq*rr or h>ui+al*qq*rr
                             for l,h,ui in zip(lo,hi,u))
    projected=max(map(abs,z))<rr and any(ui>al*rr/t0 and zi<ui-cc*rr
                                       for ui,zi in zip(u,z))
    assert original==projected,(n,rr,u,z)
    checked+=1

report={
 'kind':'exact finite algebra checks; not a replacement for general proof',
 'sympy_version':s.__version__,
 'polynomial_checks':polynomials,
 'formal_identities':{k:str(v) for k,v in identities.items()},
 'center_projection_cases':checked,
 'center_projection_parameters':{'a':'20','t':'3','dimensions':[1,2],
   'radii':['1','3/2','5/2'],'u_grid':'alpha * {0,1/4,1/2,3/4,1}',
   'z_grid':'r * {-1,-3/4,-1/2,0,1/2,3/4,1}'},
 'all_checks_passed':True,
 'r001_proof_sha256':hashlib.sha256((Path(__file__).parents[2]/'round001/theory/proof.md').read_bytes()).hexdigest()
}
assert report['r001_proof_sha256']=='6d47b46bc873021b80c9ce92fbb401b70d060626b2dcf5a28132b73e4035b551'
Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'all_checks_passed':True,'center_projection_cases':checked,
                  'polynomial_dimensions':'1..8','sympy_version':s.__version__,
                  'R001_unchanged':True}))
