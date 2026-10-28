#!/usr/bin/env python3
"""One independent exact/symbolic R090 auxiliary, not author s10 cases.

No quadrature, finite differences, grid, new tuple or dimension threshold.
The analytic proof of the local lemma and operator transfer is in review.md.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
from hashlib import sha256
import json,os,platform

HERE=Path(__file__).resolve().parent;reg=json.loads((HERE/'REGISTRATION.json').read_text())
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/audit/round090/general_s')
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')

def add(p,q):
    return [(p[i] if i<len(p) else F(0))+(q[i] if i<len(q) else F(0)) for i in range(max(len(p),len(q)))]
def scale(a,p):return [a*v for v in p]
def mul(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return r
def shift4(p):
    return [sum((p[j]*comb(j,i)*4**(j-i) for j in range(i,len(p))),F(0)) for i in range(len(p))]
def logs(y):
    z=(y-1)/(y+1);assert 0<z<1
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(16)),F(0))
    return lo,lo+2*z**33/(33*(1-z*z))

s=[F(0),F(1)];s2=mul(s,s);s4=mul(s2,s2)
Abound=scale(F(4),s);Aprime=add([F(1)],scale(F(6),s))
Dlog=add(Aprime,mul(s,Abound))
mean_x2=scale(F(4),s4);var_x2=scale(F(16),s4)
first_failure=scale(F(4),s4);source_variance=scale(F(2),s2)
Llo,Lhi=logs(F(2));Clo,_=logs(1/(1-Llo));_,Chi=logs(1/(1-Lhi))
checks={
 'new_log_derivative_polynomial_not_old_s5_constant': Dlog==[F(1),F(6),F(4)]
    and shift4(Dlog)==[F(89),F(38),F(4)] and all(v>0 for v in shift4(Dlog)),
 'new_target_moment_and_true_variance_coefficients': mean_x2==[0,0,0,0,4]
    and var_x2==[0,0,0,0,16] and first_failure==[0,0,0,0,4]
    and source_variance==[0,0,2] and F(12,11)<2,
 'uniform_local_floor_margin': F(1,2)-F(1,2)*F(4,5)==F(1,10)>0,
 'independent_limit_interval_and_sign': 0<Llo<Lhi<1
    and F(11813870,10000000)<Clo<Chi<F(11813872,10000000),
}
paths=[AUTHOR/n for n in ('PROOF.md','REGISTRATION.json','check.py','receipt.json','MANIFEST.json')]
paths += [BASE/'upper/round027/cube_cap_theory/proof.md',BASE/'revision_v3/cover/audit.md']
hashes={str(p):sha256(p.read_bytes()).hexdigest() for p in paths}
areg=json.loads((AUTHOR/'REGISTRATION.json').read_text());ar=json.loads((AUTHOR/'receipt.json').read_text())
checks['frozen_source_and_receipt_identity']=(hashes[reg['source_proof']]==reg['source_proof_sha256']
 and hashes[str(AUTHOR/'REGISTRATION.json')]==ar['registration_sha256']
 and hashes[str(AUTHOR/'check.py')]==ar['check_sha256']
 and all(sha256(Path(p).read_bytes()).hexdigest()==h for p,h in areg['sources'].items()))
assert all(checks.values()),checks
out={'status':'PASS one independent finite symbolic auxiliary; proof in review.md','pid':os.getpid(),
 'python':platform.python_version(),'author_s10_four_regimes_or_old_grid_executed':False,
 'no_numerical_integration_or_dimension_threshold':True,'checks':checks,
 'symbolic_coefficients':{'D_s':list(map(str,Dlog)),'D_s_at_s4_plus_u':list(map(str,shift4(Dlog))),
  'baseline_mean_x2_over_n':list(map(str,mean_x2)),'baseline_variance_x2_over_n':list(map(str,var_x2)),
  'first_failure_over_dmu2_n':list(map(str,first_failure)),'source_variance_x':list(map(str,source_variance))},
 'exact_intervals':{'log2_lower':str(Llo),'log2_upper':str(Lhi),'c_infty_lower':str(Clo),'c_infty_upper':str(Chi)},
 'display_only':{'c_infty_lower':float(Clo),'c_infty_upper':float(Chi),'limsup_bound_lower':float(-Chi),'limsup_bound_upper':float(-Clo)},
 'source_sha256':hashes,'registration_sha256':sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest()}
(HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','exact_intervals')},indent=2))
