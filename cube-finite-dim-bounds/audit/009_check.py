#!/usr/bin/env python3
"""Independent Q1 polynomial/jump/integral audit, SymPy exact engine."""
from pathlib import Path
from itertools import product
import hashlib,json
import sympy as s
x,y,z=s.symbols('x y z',real=True)
values=[s.Rational(11,4),s.Rational(5,8),s.Rational(1,4),s.Rational(1,8)]
poly=s.expand(sum(values[a+b+c]*(x if a else 1-x)*(y if b else 1-y)*(z if c else 1-z) for a,b,c in product((0,1),repeat=3)))
jump=s.expand(s.diff(poly,x).subs(x,0)-s.diff(poly.subs(x,-x),x).subs(x,0))
assert jump==-s.Rational(17,4)+s.Rational(7,2)*(y+z)-3*y*z
pair=s.Integer(0)
intervals=[(s.Rational(1,4),s.Rational(1,2)),(s.Rational(1,2),s.Rational(3,4))]
for (yl,yh),(zl,zh) in product(intervals,repeat=2):
 ty=(4*y-1) if yh==s.Rational(1,2) else (3-4*y)
 tz=(4*z-1) if zh==s.Rational(1,2) else (3-4*z)
 pair+=s.integrate(jump*ty*tz,(y,yl,yh),(z,zl,zh))
assert pair==-s.Rational(3,32)
mass=sum(values[sum(v!=0 for v in k)] for k in product((-1,0,1),repeat=3))
assert mass==s.Rational(21,2)
source=6*values[1]-6*values[0];assert source==-s.Rational(51,4)
out={'status':'independent exact polynomial/jump/tent integral PASS','engine':'SymPy '+s.__version__,'cell_polynomial':str(poly),'normal_jump':str(jump),'direct_piecewise_integral':str(pair),'mass':str(mass),'source_discrete_laplacian':str(source),'author_code_imported':False,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('009_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
