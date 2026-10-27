"""Exact algebra for the distribution kernel; no numerical parameter search."""
from pathlib import Path
import hashlib
import json
import sympy as s

root=Path(__file__).resolve().parent
r,q,rho,w=s.symbols('r q rho w',positive=True)
excess=s.integrate((w-q)/r,(w,q,r))
assert s.factor(excess-(r-q)**2/(2*r))==0
M=1-(r-q)**2/(4*r*rho)
assert s.factor(1-excess/(2*rho)-M)==0
# For a preface at +rho, the eligible free W interval is (-q,r).
positive_face_probability=s.integrate(1/(2*r),(w,-q,r))
assert s.factor(positive_face_probability-(r+q)/(2*r))==0
checks=['exact uniform-coordinate excess mean and normalized intersection length',
        'exact free-coordinate face capture probability']

n,mm,hh=s.symbols('n M H',positive=True)
G=mm**(n-1)/2+(n-1)*(1+q/r)*hh*mm**(n-2)/2
expanded=s.Rational(1,2)*mm**(n-1)+(n-1)*positive_face_probability*hh*mm**(n-2)
assert s.simplify(G-expanded)==0
checks.append('same-axis plus different-axis mean shell coefficient')

u,v=s.symbols('u v',real=True)
ratio=(s.exp(u)-s.exp(-v))**2/(s.exp(u)*s.exp(-v))
assert s.simplify(ratio-4*s.sinh((u+v)/2)**2)==0
checks.append('logarithmic two-point distance identity')

z,aa=s.symbols('z a',positive=True)
assert s.integrate(z*s.exp(-aa*z*z),(z,0,s.oo))==1/(2*aa)
assert s.integrate(z*s.exp(-z*z/4),(z,0,s.oo))==2
assert s.integrate(z*s.exp(-z*z/8),(z,0,s.oo))==4
checks.append('Gaussian capacity and asymptotic constants')

xx,tt=s.symbols('x t',positive=True)
assert s.factor(1/(2*((n-2)/(4*tt)))-2*tt/(n-2))==0
assert s.factor((2*tt/(n-2)).subs(tt,n/xx)-2/(xx*(1-2/n)))==0
checks.append('capacity prefactor reparameterization')
receipt={'sympy':s.__version__, 'checks':checks,
         'status':'exact finite algebra passed; the full-n theorem and limits are analytic in PROOF.md',
         'no_parameter_search':True,
         'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                           for name in ('PROOF.md','verify.py')}}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
