#!/usr/bin/env python3
"""Independent exact sign/patching checks; not a substitute for the PDE proof."""
from pathlib import Path
import hashlib,json
import sympy as s
r,z=s.symbols('r z',positive=True)
x=s.symbols('x',real=True)
rows=[]
for n in (3,4):
    kap=s.Integer(1);R=s.Integer(2);a=s.Integer(1);F=2**n
    C=kap*R**n/(n*(n-2));B=-kap*R**2/(2*n)-C*R**(2-n)
    outside=s.simplify(B+kap*r**2/(2*n)+C*r**(2-n))
    A=s.simplify(outside.subs(r,a)-(kap-F)*a**2/(2*n))
    inside=A+(kap-F)*r**2/(2*n)
    lap=lambda q:s.simplify(s.diff(q,r,2)+(n-1)*s.diff(q,r)/r)
    assert lap(inside)==kap-F and lap(outside)==kap
    assert s.simplify(inside.subs(r,a)-outside.subs(r,a))==0
    assert s.simplify(s.diff(inside,r).subs(r,a)-s.diff(outside,r).subs(r,a))==0
    assert outside.subs(r,R)==0 and s.diff(outside,r).subs(r,R)==0
    assert s.diff(inside,r).subs(r,0)==0 and A>0
    rows.append({'n':n,'f_cap':F,'kappa':1,'support_radius':2,'inner_formula':str(inside),'outer_formula':str(outside),'patching':'values and first derivatives match; radius2 derivative zero','laplacian_inside':str(lap(inside)),'laplacian_annulus':str(lap(outside))})
penalty=s.Piecewise((x*x/2,x<0),(0,True));derivative=s.diff(penalty,x)
assert derivative==s.Piecewise((x,x<0),(0,True))
assert s.simplify(z**(-s.Symbol('n'))*z**(-2)*z**s.Symbol('n')-z**(-2))==0
out={'status':'exact symbolic identities PASS','engine':'SymPy '+s.__version__,'scope':'fixed examples and sign identities only; abstract PDE proof is in008_audit.md','radial_examples':rows,'negative_part_penalty_derivative':str(derivative),'source_scaling_after_change_of_variables':'r^-2','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('008_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
