"""Predeclared finite symbolic check for next-round coarse/fine basis.
Exactly 8 open polynomial intervals + 9 endpoints, no parameter search.
Budget 30s checked per interval; no retries or alternative weights.
"""
from pathlib import Path
import sympy as s
import json,time,hashlib
r=Path(__file__).resolve().parent;out=r/'refinement_sympy.json'
if out.exists():raise RuntimeError('Do not repeat')
t=time.monotonic();x=s.Symbol('x',real=True);w=[s.Rational(k,8) for k in (1,4,6,4,1)]
def piece(y,sample):
    if sample<0:y=-y;sample=-sample
    if sample>=2:return s.Integer(0)
    if sample<=1:return s.Rational(2,3)-y**2+y**3/2
    return (2-y)**3/6
checks=[]
for a in range(-4,4):
    m=s.Rational(2*a+1,2)
    left=piece(x/2,m/2);right=sum(q*piece(x-k,m-k) for k,q in zip(range(-2,3),w))
    residual=s.Poly(s.expand(left-right),x);assert residual.is_zero
    checks.append({'interval':[a,a+1],'residual':str(residual.as_expr())});assert time.monotonic()-t<30
for a in range(-4,5):
    y=s.Integer(a);assert s.simplify(piece(y/2,y/2)-sum(q*piece(y-k,y-k) for k,q in zip(range(-2,3),w)))==0
assert sum(w)==2 and sum(w)**2==4
out.write_text(json.dumps({'status':'exact_symbolic_identity_passed','identity':'B(x/2) = sum_{k=-2}^2 w_k B(x-k)','weights':[str(q) for q in w],'polynomial_intervals':checks,'endpoint_checks':9,'outside_support':'Both sides vanish for |x|>=4','tensor_mass_factor':'4','meaning':'Exact coarse-to-fine spline conversion; does not certify any maximal-operator bound.','sympy_version':s.__version__,'seconds':time.monotonic()-t,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n');print('PASS: 8 exact polynomial pieces, 9 endpoints, tensor mass factor 4')
