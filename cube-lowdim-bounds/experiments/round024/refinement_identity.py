"""Round24 prerequisite: fixed 16 intervals + 17 endpoints, 30s guard.
Independent SymPy identity check before any LP; no search or repeat.
"""
import sympy as s
from pathlib import Path
import json,hashlib,time
root=Path(__file__).resolve().parent;out=root/'refinement_identity.json'
if out.exists():raise RuntimeError('Existing evidence; do not repeat')
t=time.monotonic();x,z=s.symbols('x z',real=True)
poly=s.Poly(s.expand((1+z+z**2+z**3)**4/64),z)
mask=[poly.nth(k) for k in range(13)]
def B(y,m):
 if m<0:y=-y;m=-m
 if m>=2:return s.Integer(0)
 if m<=1:return s.Rational(2,3)-y*y+y**3/2
 return (2-y)**3/6
res=[]
for a in range(-8,8):
 m=s.Rational(2*a+1,2)
 v=s.Poly(s.expand(B(x/4,m/4)-sum(q*B(x-k,m-k) for q,k in zip(mask,range(-6,7)))),x)
 assert v.is_zero
 res.append(str(v.as_expr()));assert time.monotonic()-t<30
for a in range(-8,9):
 a=s.Integer(a);assert s.simplify(B(a/4,a/4)-sum(q*B(a-k,a-k) for q,k in zip(mask,range(-6,7))))==0
assert sum(mask)==4
out.write_text(json.dumps({'passed':True,'mask':[str(q) for q in mask],'indices':list(range(-6,7)),'identity':'B(x/4)=sum(mask[k]*B(x-k),k=-6..6)','polynomial_residuals':res,'endpoint_checks':17,'outside':'Both sides vanish for abs(x)>=8','one_dimensional_mass_factor':4,'two_dimensional_mass_factor':16,'sympy_version':s.__version__,'seconds':time.monotonic()-t,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print('PASS: fourfold refinement, 16 pieces + 17 endpoints, exact mask',mask)
