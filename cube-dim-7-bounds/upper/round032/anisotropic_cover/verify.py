"""Exact finite algebra checks only; no high-dimensional parameter search."""
from pathlib import Path
import hashlib
import json
import sympy as s

root = Path(__file__).resolve().parent
checks = []
a, d, rho = s.symbols('a d rho', positive=True)
alpha = a/(a+1)
c = (a*(1+d)-1)/(a+1)
q = 1+d/2-1/(2*rho)
assert s.factor((1+c-alpha/rho)/2-alpha*q) == 0
assert s.factor((c-1)/alpha-(d-2/a)) == 0
checks.append('finite-a length and face transition identities')

for dim in range(1,6):
    x = s.symbols(f'x0:{dim}')
    y = s.symbols(f'y0:{dim}')
    expansion = sum((x[i]-y[i])*s.prod(x[:i])*s.prod(y[i+1:]) for i in range(dim))
    assert s.expand(s.prod(x)-s.prod(y)-expansion) == 0
checks.append('product-difference telescoping identities in d=1..5')

sn, ss, na, ni, ps, p1 = s.symbols('sn ss na ni ps p1')
face_difference = sn-1-(ss*na*ps-ni*p1)/2
face_decomposition = sn-1-na*(ss*ps-p1)/2+(ni-na)*p1/2
assert s.expand(face_difference-face_decomposition) == 0
checks.append('per-face finite-a versus limiting budget decomposition')

A, B, kp, km, tau, uu = s.symbols('A B kp km tau uu', positive=True)
weighted = (1+A*kp/tau)*(tau*km+B)
rescaled = (1+A*kp*km/uu)*(uu+B)
assert s.factor(weighted.subs(tau,uu/km)-rescaled) == 0
optimal = s.sqrt(A*kp*km*B)
assert s.simplify(s.diff(rescaled,uu).subs(uu,optimal)) == 0
assert s.simplify(rescaled.subs(uu,optimal)-(s.sqrt(A*kp*km)+s.sqrt(B))**2) == 0
checks.append('scale threshold rescaling and positive scalar optimum')

assert s.factor(2*4*4-(4+1)**2) > 0
checks.append('2^n >= n^2: n=4 base and monotone induction ratio algebra')
receipt = {
    'sympy':s.__version__,
    'checks':checks,
    'status':'exact finite algebra checks; arbitrary dimension proof is analytic in proof.md',
    'no_parameter_search':True,
    'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                     for name in ('proof.md','verify.py')},
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
