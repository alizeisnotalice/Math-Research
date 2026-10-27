#!/usr/bin/env python3
"""Independent symbolic check of the primitive and finite-a identities."""
from pathlib import Path
import datetime, hashlib, json, sys, time
import sympy as s

D=Path(__file__).resolve().parent
started=datetime.datetime.now(datetime.timezone.utc).isoformat()
wall=time.perf_counter()
x,p,t,a,rho=s.symbols('x p t a rho',positive=True)
polynomial=s.Poly(s.expand((p-x/2)**7),x)
primitive=p**7*s.log(x)+sum(polynomial.coeff_monomial(x**j)*x**j/j for j in range(1,8))
assert s.simplify(s.diff(primitive,x)-(p-x/2)**7/x)==0
alpha=a/(a+1)
c=(a*(1+1/t)-1)/(a+1)
expressions={
    'c_plus_one':1+c-alpha*(2+1/t),
    'c_minus_one':c-1-(a/t-2)/(a+1),
    'finite_corner_kernel':(1+c-alpha/rho)/2-alpha*(1+1/(2*t)-1/(2*rho)),
    'corner_inner':(1+c-alpha)/2-alpha*(1+1/t)/2,
    'corner_transition_gap':alpha/(c-1)-t-2*t*t/(a-2*t),
}
for name,value in expressions.items():
    assert s.factor(value)==0,name
integrated=[s.factor(polynomial.coeff_monomial(x**j)/p**(7-j)/j) for j in range(1,8)]
assert integrated==list(map(s.Rational,['-7/2','21/8','-35/24','35/64','-21/160','7/384','-1/896']))
result={
    'status':'PASS','primitive_derivative_zero':True,
    'polynomial':str(polynomial.as_expr()),
    'integrated_coefficients':list(map(str,integrated)),
    'five_finite_a_identities':{name:str(s.factor(value)) for name,value in expressions.items()},
    'start_utc':started,'finish_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'wall_seconds':time.perf_counter()-wall,'python':sys.version,'executable':sys.executable,
    'sympy_version':s.__version__,
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(D/'symbolic_result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
