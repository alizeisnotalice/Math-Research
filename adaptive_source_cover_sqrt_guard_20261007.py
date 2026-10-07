"""Exact guard for eta0/sqrt(n) traffic policy; no old data rerun."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json
import hashlib

out=Path(__file__).with_suffix('.json')
assert not out.exists(),'Refuse overwrite.'
eta0=F(49,65536);epsilon=F(49,16384)
assert 4*eta0==epsilon
rows=[]
for n in (4,16,64,8,32,128):
    d=2**20
    p=isqrt(n*d*d)
    lo=F(p,d)
    exact=(p*p==n*d*d)
    hi=lo if exact else F(p+1,d)
    assert lo*lo<=n<=hi*hi and hi-lo<=F(1,d)
    eta_safe=eta0/hi
    eta_hi=eta0/lo
    assert 4*eta_safe<=epsilon
    paid_prefactor_without_4e6=1/eta_safe
    assert paid_prefactor_without_4e6==hi/eta0
    assert eta_safe<=eta_hi and hi<=lo+F(1,d)
    rows.append({'n':n,'sqrt_lower':str(lo),'sqrt_upper':str(hi),'sqrt_exact_rational':exact,
                 'eta_conservative_rational':str(eta_safe),'eta_irrational_policy_upper_bound':str(eta_hi),
                 'near_absorption_coefficient_exact':str(4*eta_safe),'epsilon':str(epsilon),
                 'paid_traffic_factor_without_4e6_exact':str(paid_prefactor_without_4e6),
                 'max_paid_extra_factor_additive_without_4e6':str(F(1,d)/eta0)})
res={'status':'PASS_EXACT_SQRT_TRAFFIC_POLICY_GUARD','eta0':str(eta0),'rows':rows,
     'policy':'h=2a/n, eta_n=eta0/sqrt n; paid traffic<=4e6 sqrt n W/eta0. Conservative rational guards round eta down, not up.',
     'remaining':'Same-receiver far traffic unproved; no general full-energy closure.',
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out.write_text(json.dumps(res,indent=2))
print(json.dumps({'status':res['status'],'dimensions':[r['n'] for r in rows]}))
