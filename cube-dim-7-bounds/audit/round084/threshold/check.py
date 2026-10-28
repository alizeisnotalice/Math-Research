#!/usr/bin/env python3
"""Independent R084 N12/T8 audit; no execution/import of producer 55 gates."""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import hashlib,json,os,platform

HERE=Path(__file__).resolve().parent
SOURCE=Path('/tmp/r079-compute-https-20261028/cube-dim-7-bounds/upper/round084/threshold_compute')
PROOF=Path('/Users/zhengzhihao/Desktop/T/cube-r078-theory-bernstein-sparse-20261028/cube-dim-7-bounds/upper/round083/general_s/PROOF.md')
own=json.loads((HERE/'REGISTRATION.json').read_text())
L=F(own['L0']);sqrtL=F(2**25);N=L**12/F(factorial(12));s=F(9,2)
x=F(387,1000)*L;h=1/L**2;e=F(1,100000);sig=F(613,1000);chi=F(1583,2000)
r0=1+F(1,1000000);r1=F(4499549,1000000);q1=F(9999,10000)
rho=N/(s*x);mean=49*s*s*x*x/N+14*s*x/N;var=196*s*s*x*x/N
P0=1-2*e;ulo=P0*(1-e)/(1+e);uhi=(1+e)/(P0*(1-e));M=2*(1+e)*uhi
b=2*F(99999,100000)*(1-e)*ulo-2*e;Gamma=F(1503,1000)*b/M
nu=1/s-(1+1/s)*e;eta=(1-q1)/(2*q1);K=(nu-1/N)*eta*x/(1+e)
gap=F(3,100000)*x-F(3,5);bad=2*x/gap**2;cT=sig/2-chi*(1-sig)
T8=lambda y:y**8/F(factorial(8))
cap=M/T8(K);tau=r1*L**2*sqrtL/T8(cT*L)
root=1/N+1/T8(sig*L)+2*sqrtL/T8(sig*L/2)

def logs(v):
    z=(v-1)/(v+1)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)),F(0))
    return lo,lo+2*z**41/(41*(1-z*z))

masslo,_=logs(r1/r0);_,pricehi=logs(F(1000,387));Qlo,Qhi=logs(F(18,11))
muhi=F(7,2)-F(11,2)*Qlo
true_y=[-F(5,2),-F(3,2),-F(1,4),F(1,12),F(1,16)]
weak_y=[-F(5,2),-F(3,2),-F(1,2),F(1,12),F(1,16)]
shift=lambda cc:[sum((cc[j]*comb(j,i)*4**(j-i) for j in range(i,5)),F(0)) for i in range(5)]
true_u=shift(true_y);weak_u=shift(weak_y)
checks={
 'fixed_base_explicit_geometry_and_local_domain': L==2**50 and sqrtL**2==L and N>=506250
    and x>=2 and x/N<1 and x*N*N>1 and 2*h<F(1,1000000)
    and h<=F(1,1000) and r0-2*h>1 and 1-h>0 and (s-1)/s-h>0
    and 1/(x*N*N)<1 and q1/(x*N*N)<=h/s and N/x>=50
    and N/x>r1 and rho>4*r1 and rho>=14 and r1<5 and 4*r1>6 and 1+r1<6,
 'new_baseline_count_target_cap_and_profile_probabilities': mean<=e/2 and 4*var/e**2<=e/8
    and 4/(e*e*N)<=e/8 and 1+x/(2*N)-1/(2*r1)<=F(9,10)
    and F(9,10)**299<=e/8 and F(81,16)/N<=F(1,100000)
    and F(64072265625)/N<=e/8,
 'whole_matching_new_roundings_and_narrow_retention_gap': 7/rho<=e and 7/rho<=F(1,2)
    and 180*x*h<=e and 2*h+1/N<=3*e/4 and 80*h<=F(1,100000)
    and F(79144,100000)/(1-e)<F(79147,100000) and F(27,16)/(1-e)<2
    and gap>0 and bad<=e,
 'count_two_ninths_weighted_cap_and_tau_T8': nu==F(22221,100000) and eta==F(1,19998)
    and nu-1/N>0 and K>=4 and cap<e and 1/(1-(h/r0)**2)<=1+e
    and cT==F(379,2000000) and cT*L>=4 and tau<e,
 'new_fixed_mass_deficit_profile_price': Qhi<F(1,2) and muhi<F(7914,10000)
    and (1-e)*masslo>F(1503,1000) and (1-e)*Gamma>F(1501,1000)
    and pricehi<F(19,20) and F(1501,1000)-F(19,20)-F(11,20)==F(1,1000),
 'positive_root_complete_cost_and_actual_n_conditions': L**2/8>=x and L/4>F(1501,1000)
    and sig*L/2>=2 and root<F(1,1000) and N>=2 and rho>6,
 'both_true_and_weaker_quartic_tail_polynomials': true_u==[F(53,6),F(33,2),F(27,4),F(13,12),F(1,16)]
    and weak_u==[F(29,6),F(29,2),F(13,2),F(13,12),F(1,16)]
    and all(v>0 for v in true_u+weak_u),
}
reg=json.loads((SOURCE/'registration.json').read_text());receipt=json.loads((SOURCE/'result.json').read_text())
paths=[SOURCE/n for n in ('registration.json','README.md','verify.py','result.json')]+[PROOF]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
checks['source_registration_script_and_frozen_proof_identity']=(F(reg['L0'])==L
 and hashes[str(PROOF)]==reg['proof_sha256']==receipt['proof_sha256']==own['proof_sha256']
 and hashes[str(SOURCE/'registration.json')]==receipt['registration_sha256']
 and hashes[str(SOURCE/'verify.py')]==receipt['script_sha256'])
assert all(checks.values()),checks
values={'L0':L,'N12':N,'rho_lower':rho,'baseline_mean_upper':mean,
 'baseline_failure_upper':4*var/e**2,'count_failure_upper':4/(e*e*N),
 'first_profile_failure_upper':F(64072265625)/N,'chebyshev_gap':gap,'face_bad_upper':bad,
 'weighted_cap_exponent_lower':K,'cap_T8_upper':cap,'tau_T8_upper':tau,'root_T8_upper':root,
 'mass_margin':(1-e)*masslo-F(1503,1000),'paid_gamma_margin':(1-e)*Gamma-F(1501,1000)}
out={'status':'PASS threshold conditions; conditional on frozen R083 full-chain acceptance',
 'pid':os.getpid(),'python':platform.python_version(),'original_55_executed':False,
 'recipe':'single registered L0; Fraction, N12, positive T8, log20 and exact quartic shifts',
 'checks':checks,'source_sha256':hashes,'own_registration_sha256':hashlib.sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest(),
 'exact_bounds':{k:str(v) for k,v in values.items()},'display_only':{k:float(v) for k,v in values.items()},
 'true_quartic_shift':[str(v) for v in true_u],'weaker_quartic_shift':[str(v) for v in weak_u]}
(HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','exact_bounds')},indent=2))
