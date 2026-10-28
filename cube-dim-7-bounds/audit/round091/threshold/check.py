#!/usr/bin/env python3
"""Independent R091 fixed-threshold arithmetic, not author58 or worker15.

Weaker N12/T8 and binary range-reduced log20 enclosures differ from
the producer N16/S24/log96 recipe. One fixed threshold only.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from hashlib import sha256
import json,os,platform

HERE=Path(__file__).resolve().parent;own=json.loads((HERE/'REGISTRATION.json').read_text())
SOURCE=Path('/tmp/r079-compute-https-20261028/cube-dim-7-bounds/upper/round091/general_s_model')
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
R090=BASE/'upper/round090';V={k:F(v) for k,v in own['parameters'].items()}

def rawlog(y):
    if y==1:return F(0),F(0)
    z=(y-1)/(y+1);assert 0<z<F(1,3)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)),F(0))
    return lo,lo+2*z**41/(41*(1-z*z))

# Endpoint log2 is computed directly, then every other argument lies in [1,2).
z=F(1,3);L2lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)),F(0))
L2hi=L2lo+2*z**41/(41*(1-z*z))
def logs(y):
    assert y>1
    k=0
    while y>=2:y/=2;k+=1
    lo,hi=rawlog(y)
    return k*L2lo+lo,k*L2hi+hi

L=F(own['L0']);N=L**12/F(factorial(12));s=V['s'];chi=V['chi'];sig=V['sigma'];d=V['d_tau']
x=(1-sig)*L;h=1/L**2;e=V['e'];r0=1+V['trim'];r1=s*V['qmax']-V['trim'];rho=N/(s*x)
Qlo,Qhi=logs(2*s/(s+1));mulo=s-1-(s+1)*Qhi;muhi=s-1-(s+1)*Qlo
dlo=(chi-muhi)/4;dhi=(chi-mulo)/4;meanface_hi=(muhi+chi)/(2*(1-e))
masslo,_=logs(r1/r0);rawlo,_=logs(s*V['qmax']);_,pricehi=logs(1/(1-sig))
P0=1-2*e;ulo=P0*(1-e)/(1+e);uhi=1/ulo;M=2*(1+e)*uhi
b=2*V['Ret']*(1-e)*ulo-2*e;Gamma=V['m']*b/M;slack=V['gamma']-V['price']-V['c']
Ds=4*s*s+6*s+1;mean=4*s**4*x*x/N+4*s*s*x/N;var=16*s**4*x*x/N
nu=1/s-(1+1/s)*e;eta=(1-V['qmax'])/(2*V['qmax']);K=(nu-1/N)*eta*x/(1+e)
aface=chi-V['m_face'];gap=aface*x-1;facebad=2*s*s*x/gap**2
ct=d-chi*(1-sig);crossrate=sig-d
T8=lambda y:y**8/F(factorial(8))
targetcap=1/T8(F(10000)/(4*s));cap=M/T8(K);tau=r1*L**2/T8(ct*L)
Hrem=1/N+1/T8(d*L)+1/T8(sig*L)+L/T8(crossrate*L)
checks={
 'fixed_tuple_s_dependent_domains_and_scalar_gaps': s>=4 and 0<sig<1 and 1/s<V['qmax']<1
    and 0<V['gamma']<rawlo and V['c']<V['gamma']-pricehi and 0<mulo<muhi<chi
    and chi*(1-sig)<d<sig and 0<e<min(F(1,12),1/(2*(s+1)))
    and (1-e)*masslo>V['m']>V['gamma'] and V['Ret']*V['m']>V['gamma']
    and meanface_hi<V['m_face']<chi and (1-e)*Gamma>V['gamma']
    and pricehi<V['price'] and slack==F(79131,100000000)>0,
 'uniform_arbitrary_s_local_and_complete_matching_domain': L==2**33 and N>=10001
    and 0<h<=F(1,1000) and 2*h<V['trim'] and r0<r1 and r0-2*h>1
    and 1-h>0 and (s-1)/s-h>0 and 1/(x*N*N)<min(1,h)
    and V['qmax']/(x*N*N)<=h/s and x>=2 and x/N<1 and x*N*N>1
    and N/x>=2*r1 and rho>4*(s+1) and 1+r1<s+1
    and L2hi<1 and 2*s/rho<F(1,2) and 2*s/rho<=e,
 'five_target_events_new_s4_moments_and_cap_T8': mean<=e/2 and 4*var/e**2<=e/8
    and 4/(e*e*N)<=e/8 and 1+x/(2*N)-1/(2*r1)<=1-1/(4*s)
    and targetcap<e/8 and dlo>0 and s*s/N<=dlo/2
    and 4*s**4/(dlo*dlo*N)<=e/8,
 'whole_matching_Ds_5s_and_true_source_retention': 2*Ds*x*h<=e and 10*s*h<=dlo
    and 2*h+1/N<=3*e/4 and s*s/(1-e)<2*s*s
    and aface>0 and gap>0 and facebad<=1-V['Ret'],
 'weighted_cap_and_free_tau_complete_payments_T8': nu>1/(2*s) and nu-1/N>0 and eta>0
    and K>=4 and cap<e and 1/(1-(h/r0)**2)<=1+e
    and ct>0 and ct*L>=4 and tau<e,
 'complete_H_not_F_and_actual_history_size_conditions': L**2/8>=x and L/4>V['gamma']
    and min(sig,d,crossrate)*L>=2 and Hrem<slack and N>=2 and rho>s+1,
}
paths=[SOURCE/n for n in ('exact_report.md','exact_verify.py','exact_registration.json',
 'exact_result_L8589934592.json','exact_manifest.json','summary.json')]
paths += [R090/'general_s/PROOF.md',R090/'PRIMARY_ACCEPTANCE.md',R090/'REPORT.md',
 BASE/'audit/round090/full/review.md']
hashes={str(p):sha256(p.read_bytes()).hexdigest() for p in paths}
reg=json.loads((SOURCE/'exact_registration.json').read_text());res=json.loads((SOURCE/'exact_result_L8589934592.json').read_text())
checks['candidate_and_accepted_R090_hash_identity']=(hashes[str(SOURCE/'exact_registration.json')]==own['source_exact_registration_sha256']==res['registration_sha256']
 and hashes[str(SOURCE/'exact_verify.py')]==res['script_sha256']
 and hashes[str(SOURCE/'summary.json')]==reg['source_summary_sha256']
 and hashes[str(R090/'general_s/PROOF.md')]==reg['proof_sha256']==own['R090_proof_sha256'])
assert all(checks.values()),checks
vals={'N12':N,'s':s,'rho_lower':rho,'D_s':Ds,'mu_lower':mulo,'mu_upper':muhi,'dmu_lower':dlo,'dmu_upper':dhi,
 'source_mean_upper':meanface_hi,'log_transport_upper':2*Ds*x*h,'target_cap_T8_upper':targetcap,
 'baseline_failure_upper':4*var/e**2,'count_failure_upper':4/(e*e*N),
 'profile_failure_upper':4*s**4/(dlo*dlo*N),'source_face_bad_upper':facebad,'allowed_source_bad':1-V['Ret'],
 'weighted_exponent_lower':K,'weighted_cap_T8_upper':cap,'tau_rate':ct,'H_cross_rate':crossrate,
 'tau_T8_upper':tau,'full_H_remainder_T8_upper':Hrem,'mass_margin':(1-e)*masslo-V['m'],
 'paid_gamma_margin':(1-e)*Gamma-V['gamma'],'scalar_slack':slack}
out={'status':'PASS independent exact threshold arithmetic; full-tail/operator acceptance in review.md',
 'pid':os.getpid(),'python':platform.python_version(),'recipe':'single L0=2^33; N12/T8 and binary range-reduced log20',
 'author58_or_worker15_executed':False,'other_threshold_evaluated':False,'checks':checks,
 'exact_bounds':{k:str(v) for k,v in vals.items()},'display_only':{k:float(v) for k,v in vals.items()},
 'source_sha256':hashes,'registration_sha256':sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest()}
(HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('exact_bounds','source_sha256')},indent=2))
