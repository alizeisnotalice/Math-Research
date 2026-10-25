#!/usr/bin/env python3
"""Exact finite checks for analytic lemmas; not a C5 bound certificate."""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from flint import arb,ctx
start=time.monotonic()
n=5;q=s.Rational(5,4)
B=s.ones(5)-2*s.eye(5);B[:,0]=-B[:,0]
assert abs(B.det())==48 and all(B[0,j]==1 for j in range(5))
V=s.Matrix([[1,*v] for v in itertools.product([-1,1],repeat=4)])
assert V.T*V==16*s.eye(5)
rate=sum(abs(sum(e)) for e in itertools.product([-1,1],repeat=4))
assert rate==24
exponent=s.Rational(-3,5)+s.Rational(4,5)
assert exponent==s.Rational(1,5)
c,g=s.symbols('c g',positive=True)
f=c**(-s.Rational(1,4))/(g-c)
assert s.simplify(s.diff(f,c)/f-(5*c-g)/(4*c*(g-c)))==0
ctx.prec=256
Cf=arb(256)/(3*arb.pi()**5)
A=arb(16)
Ecoeff=(arb(5)/6)*Cf.root(5)*(A/arb.pi()).root(5)**4/(arb(16).root(5)**4)
records=[]
for k in [4,6,8,10]:
 eps=arb(1)/2**k
 E=Ecoeff*eps.root(5)**6
 records.append({'epsilon':f'1/{2**k}','norm_budget_ball':str(E)})
res={'status':'finite_analytic_checks_passed_not_new_C5_bound','q':'5/4','basis_det_abs':'48','Gram':'16 I_5','free_coordinate_symbol_lower_factor':'6','face_sup_coefficient':'256/(3*pi^5)','whole_face_marginal':'16*t/(pi*(1+256*t^2))','short_time_exponent':'1/5','short_time_integral_exponent':'6/5','positive_face_short_time_norm_coefficient':'(5/6)*(256/3)^(1/5)*pi^(-9/5)','scalar_optimum':'c=g/5','cost_fourth_power':'3125*I0^4/g^5','illustrative_outward_norm_balls':records,'seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy':s.__version__,'arb_precision':ctx.prec}
print(json.dumps(res,ensure_ascii=False,indent=2))
