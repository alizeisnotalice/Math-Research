"""Exact algebra checks for the compact radial supersolution, not a PDE solver."""
import hashlib,json
from pathlib import Path
import sympy as s

r,a,R,k=s.symbols('r a R k',positive=True)
rows=[]
for n in (3,4):
    F=k*(R/a)**n
    outer=F*a**n/(n*(2-n))*(R**(2-n)-r**(2-n))-k/(2*n)*(R**2-r**2)
    inner=outer.subs(r,a)+(F-k)/(2*n)*(a*a-r*r)
    lap=lambda u:s.simplify(s.diff(u,r,2)+(n-1)/r*s.diff(u,r))
    checks=[s.simplify(-lap(inner)-(F-k)),s.simplify(-lap(outer)+k),
            s.simplify((inner-outer).subs(r,a)),
            s.simplify((s.diff(inner,r)-s.diff(outer,r)).subs(r,a)),
            s.simplify(outer.subs(r,R)),s.simplify(s.diff(outer,r).subs(r,R))]
    assert all(z==0 for z in checks)
    rows.append({'n':n,'six_exact_residuals':list(map(str,checks))})
print(json.dumps({'rows':rows,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
