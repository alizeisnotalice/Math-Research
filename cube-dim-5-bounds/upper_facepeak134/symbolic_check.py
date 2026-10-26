"""Registered finite NEW near/far geometry identities only; no C endpoints."""
from pathlib import Path
import hashlib
import json
import os
import time
import sympy as s

HERE = Path(__file__).resolve().parent
assert not (HERE/'symbolic_results.json').exists()
regraw = (HERE/'symbolic_registration.json').read_bytes()
reg = json.loads(regraw)
assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == reg['code_sha256']
for row in reg['source_locks']:
    assert hashlib.sha256(Path(row['path']).read_bytes()).hexdigest() == row['sha256']
started = time.monotonic()
r,m = s.symbols('r m',real=True)
f = 5-4*s.sqrt(5)*r+5*r*r
assert s.expand(5*(1-2*r/s.sqrt(5))**2+r*r-f) == 0
assert s.simplify(s.diff(f,r).subs(r,2/s.sqrt(5))) == 0
assert s.diff(f,r,2) == 10
assert s.Rational(1,10)**2 < s.Rational(4,5)
assert 5-9*s.Rational(1,10)+5*s.Rational(1,10)**2 == s.Rational(83,20)
assert s.Rational(16,480) == s.Rational(1,30)
assert s.Rational(9,4)**2 > 5
P = s.eye(5)-s.ones(5)/5
assert (P*P-P) == s.zeros(5)
assert (P[:,0].dot(P[:,0])) == s.Rational(4,5)
# General m KKT residuals; no old source geometry enumeration or baseline.
positive_count = 4-m
h = 1/(m+1)
b = s.cancel(1+positive_count*h)
Qstar = 1+positive_count*h*h-b*b/5
assert s.cancel(Qstar-m/(m+1)) == 0
assert s.cancel(2*h-2*b/5) == 0
assert s.cancel(2*b/5-2/(m+1)) == 0
result = {'status':'new_near_far_symbolic_identities_passed','pid':os.getpid(),
    'seconds':time.monotonic()-started,'sympy':s.__version__,
    'registration_sha256':hashlib.sha256(regraw).hexdigest(),
    'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks':['near_norm_expansion','near_polynomial_monotonicity_endpoint',
      'eta1over10_domain_and_d83over20','16over480_equals1over30',
      'sqrt5_less9over4','projection_e1_norm_squared4over5',
      'general_m_face_KKT_residuals'],
    'real_C_endpoint_arithmetic':0,'new_H_gradient_log_Gamma_Arb_calls':0,
    'old_baseline_runs':0,'matlab_calls':0}
with (HERE/'symbolic_results.json').open('x') as f:
    json.dump(result,f,ensure_ascii=False,sort_keys=True,indent=2); f.write('\n')
print(json.dumps(result,ensure_ascii=False))
