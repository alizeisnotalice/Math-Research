"""Exact formula and bounded numerical checks for L_V (1-|x|^2)_+^2."""
from pathlib import Path
import json,hashlib,math,itertools
import sympy as sy
import numpy as np
from scipy.integrate import quad
s,a,h,n=sy.symbols('s a h n',real=True)
f=n*n*(h*h-(s-a)**2)**2
q,r=sy.div(sy.diff(f,s),s,s)
assert sy.simplify(q-4*n*n*(s*s-3*a*s+3*a*a-h*h))==0
assert sy.simplify(r-4*n*n*a*(h*h-a*a))==0
poly=sy.integrate(q,(s,a-h,a+h));assert sy.simplify(poly-4*n*n*(2*a*a*h-sy.Rational(4,3)*h**3))==0

def one(x,v):
 x=np.array(x,dtype=float);v=np.array(v,dtype=float);n=float(v@v);a=-float(x@v)/n;A=1-float(x@x)+float(x@v)**2/n
 if A<=0:return 0.,0.,'line misses or is tangent to ball'
 h=math.sqrt(A/n);a2=a*a;h2=h*h
 if a==0:logterm=0.
 elif abs(abs(a)-h)<1e-13:logterm=0.
 else:logterm=a*(h2-a2)*math.log(abs((a+h)/(a-h)))
 answer=4*n*n/math.pi*(2*a2*h-4*h**3/3+logterm)
 # PV integration of f'(s)/s by QUADPACK cauchy rule, with endpoint limiting case exact polynomial remainder.
 def fp(s):return -4*n*n*(s-a)*(h*h-(s-a)**2)
 if a-h<0<a+h:
  numeric=quad(fp,a-h,a+h,weight='cauchy',wvar=0,epsabs=1e-11)[0]/math.pi
 elif abs(abs(a)-h)<1e-13:
  numeric=quad(lambda s:4*n*n*(s*s-3*a*s+3*a*a-h*h),a-h,a+h,epsabs=1e-11)[0]/math.pi
 else:numeric=quad(lambda s:fp(s)/s,a-h,a+h,epsabs=1e-11)[0]/math.pi
 return answer,numeric,'PV derivative integral'
rows=[]
for dim in (3,4):
 pts=[[0]*dim,[.2]+[.1]*(dim-1),[1]+[0]*(dim-1),[2]*dim,[3,0]+[0]*(dim-2)]
 vv=[(1,)+v for v in itertools.product((-1,1),repeat=dim-1)]
 for point in pts:
  vals=[one(point,v) for v in vv];err=max(abs(a-b) for a,b,_ in vals);assert err<1e-10
  rows.append({'n':dim,'x':point,'Lb_formula':sum(a for a,b,_ in vals),'Lb_quad':sum(b for a,b,_ in vals),'max_direction_error':err})
 radius=sy.symbols('radius',nonnegative=True)
 mass=sy.simplify(2*sy.pi**(sy.Rational(dim,2))/sy.gamma(sy.Rational(dim,2))*sy.integrate((1-radius**2)**2*radius**(dim-1),(radius,0,1)))
 rows.append({'n':dim,'b_mass':str(mass)})
out={'status':'formula algebra and finite quadrature checks PASS; no corrected kernel certified','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'scope':'Lb is full directional generator including tails; finite checks not full source dominance','shared_libraries':['SymPy','SciPy','numpy']}
(Path(__file__).resolve().parents[1]/'results/round009_bump_interface.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
