#!/usr/bin/env python3
"""New targeted R080 audit, not a replay of the producer 53/author 18 gates.

N16, log20 and P40 independently check the new weighted deletion margin,
the selected gamma/scalar slack, and exact complete-old-cost separators.
All source inputs are read-only. No producer program is imported/executed.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,os,platform

HERE=Path(__file__).resolve().parent
THEORY=Path('/Users/zhengzhihao/Desktop/T/cube-r078-theory-bernstein-sparse-20261028/cube-dim-7-bounds/upper')
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/upper')

def log20(y):
    z=(y-1)/(y+1)
    assert 0<z<1
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)),F(0))
    return lo,lo+2*z**41/(41*(1-z*z))

def P40(y):
    assert y>0
    return sum((y**j/F(factorial(j)) for j in range(41)),F(0))

L=F(4096);N=L**16/F(factorial(16));x=F(89,100)*L
q=F(247,250);delta=F(1,2**19);eJ=eR=eS=eD=F(1,2**20)
eN=eB=eT=delta;epsA=F(1,2**16);D=F(1536)
P0=1-eB-6/D;ulo=P0*(1-eN)/(1+eJ);uhi=(1+eN)/(P0*(1-eJ))
M=2*(1+eR)*uhi;nu=F(2,3)-F(5,3)*eN;eta=(1-q)/(2*q)
K=(nu-1/N)*eta*x/(1+eJ);cap=M/P40(K)
rhi=F(3,2)*q-delta;rlo=1+delta
ml=F(196691,500000);gamma=F(39029,100000);price=F(5827,50000);c=F(5473,20000)
ret=F(999999,1000000);b=2*ret*(1-eR)*ulo-eS-eT;Gamma=ml*b/M
masslo,_=log20(rhi/rlo);_,pricehi=log20(F(100,89))
scalar=1/N+1/P40(F(11,100)*L)+128/P40(F(11,200)*L)
margin69=c-F(234,1000)-F(4,2**30)
margin78=c-F(599,2000)+F(25,256)-F(4,8192)

checks={
 'new_weighted_deletion_with_arbitrary_pin': nu>0 and nu-1/N>0 and eta==F(3,494)
    and K>4 and cap<eS,
 'new_mass_gamma_and_complete_scalar': (1-epsA)*masslo>ml and 0<gamma<(1-eD)*Gamma<1
    and pricehi<price and gamma-price-c==F(1,10000) and scalar<F(1,10000),
 'complete_old_cost_separators': margin69>0 and margin78>0,
 'obsolete_unweighted_rate_is_negative': (eta-F(51,1000))*F(89,100)<0,
}
files=[THEORY/'round080/weighted_threshold/PROOF.md',THEORY/'round080/weighted_threshold/verify.py',
 THEORY/'round080/weighted_threshold/verification.json',THEORY/'round080/weighted_threshold/REGISTRATION.json',
 THEORY/'round079/weighted_cap/PROOF.md',MAIN/'round080/compute/case_988_1000.json',
 MAIN/'round069/profile_threshold/THRESHOLD_PROOF.md']
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
reg=json.loads(files[3].read_text());receipt=json.loads(files[2].read_text())
checks['read_only_source_identity']=hashes[str(files[5])]==reg['producer_locks']['case_988_1000.json']
checks['compact_receipt_registration_identity']=(hashes[str(files[3])]==receipt['registration_sha256'])
assert all(checks.values()),checks
values={'N16':N,'nu':nu,'K_lower':K,'weighted_cap_upper':cap,'weighted_cap_slack':eS-cap,
 'paid_gamma_margin':(1-eD)*Gamma-gamma,'new_mass_margin':(1-epsA)*masslo-ml,
 'scalar_remainder_upper':scalar,'price_margin':price-pricehi,'Phi69_gap_lower_coefficient':margin69,
 'Phi78_gap_lower_coefficient':margin78}
out={'status':'PASS targeted sensitive arithmetic; proof audit in review.md','pid':os.getpid(),
 'python':platform.python_version(),'method':'exact Fraction; N16, log20 positive atanh+tail, P40',
 'producer_53_or_compact_18_executed':False,'checks':checks,'source_sha256':hashes,
 'exact':{k:str(v) for k,v in values.items()},'display_only':{k:float(v) for k,v in values.items()}}
(HERE/'sensitive_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'checks':checks,'display_only':out['display_only']},indent=2))
