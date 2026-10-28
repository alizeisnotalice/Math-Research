#!/usr/bin/env python3
"""Independent R087 fixed-tuple audit; no import/execution of author 8+56 gates."""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import hashlib,json,os,platform

HERE=Path(__file__).resolve().parent
SOURCE=Path('/tmp/r079-compute-https-20261028/cube-dim-7-bounds/upper/round087/threshold_candidate')
PROOF=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/upper/round085/parameter_family/PROOF.md')
AUDIT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/audit/round085/full')
reg=json.loads((HERE/'REGISTRATION.json').read_text());P={k:F(reg[k]) for k in
 ('L0','s','chi','sigma','qmax','gamma','c','trim','e','Ret','m','m_face','pin','price','slack')}
L=P['L0'];s=P['s'];chi=P['chi'];sig=P['sigma'];q=P['qmax'];gamma=P['gamma'];c=P['c']
e=P['e'];ret=P['Ret'];mass=P['m'];face=P['m_face'];pin=P['pin'];price=P['price'];slack=P['slack']
sqrtL=F(2**40);N=L**12/F(factorial(12));x=(1-sig)*L;h=1/L**2
r0=1+P['trim'];r1=s*q-P['trim'];rho=N/(s*x)

def logs(v):
    z=(v-1)/(v+1)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(32)),F(0))
    return lo,lo+2*z**65/(65*(1-z*z))

Qlo,Qhi=logs(2*s/(s+1));mulo=s-1-(s+1)*Qhi;muhi=s-1-(s+1)*Qlo
dlo=(chi-muhi)/4;dhi=(chi-mulo)/4
# Actual mu+2d=(mu+chi)/2, so this upper bound preserves correlation.
transport_hi=(muhi+chi)/2
masslo,_=logs(r1/r0);bandlo,_=logs(s*q);_,pricehi=logs(1/(1-sig))
P0=1-2*e;ulo=P0*(1-e)/(1+e);uhi=(1+e)/(P0*(1-e));M=2*(1+e)*uhi
b=2*ret*(1-e)*ulo-2*e;Gamma=mass*b/M;ct=sig/2-chi*(1-sig)
mean=81*s*s*x*x/N+18*s*x/N;var=324*s*s*x*x/N
nu=1/s-(1+1/s)*e;eta=(1-q)/(2*q);K=(nu-1/N)*eta*x/(1+e)
gap=(chi-face)*x-pin;bad=10*x/gap**2
T8=lambda y:y**8/F(factorial(8))
cap=M/T8(K);tau=r1*L**2*sqrtL/T8(ct*L)
root=1/N+1/T8(sig*L)+2*sqrtL/T8(sig*L/2)
pc=[-F(5,2),-F(3,2),-F(1,2),F(1,12),F(1,16)]
shift=[sum((pc[j]*comb(j,i)*4**(j-i) for j in range(i,5)),F(0)) for i in range(5)]
checks={
 'strict_primary_parameters_and_secondary_choice': 4<=s<=5 and 0<sig<1 and 1/s<q<1
    and chi>muhi and ct>0 and 0<gamma<bandlo and c<gamma-pricehi
    and 0<e<F(1,12) and (1-e)*masslo>mass>gamma and ret*mass>gamma
    and transport_hi/(1-e)<face<chi and (1-e)*Gamma>gamma
    and pricehi<price and gamma-price-c==slack>0,
 'explicit_uniform_local_and_full_matching_domain': L==2**80 and sqrtL**2==L and N>=400
    and x>=2 and x/N<1 and x*N*N>1 and 2*h<P['trim'] and h<=F(1,1000)
    and r0-2*h>1 and 1-h>0 and (s-1)/s-h>0 and 1/(x*N*N)<1
    and q/(x*N*N)<=h/s and N/x>=2*r1 and rho>20 and rho>4*r1
    and 1<r1<5 and Qhi<pin and 9/rho<=F(1,2) and 9/rho<=e,
 'five_event_probability_budgets_with_actual_dmu_lower': dlo>0 and dlo<=dhi
    and mean<=e/2 and 4*var/e**2<=e/8 and 4/(e*e*N)<=e/8
    and 1+x/(2*N)-1/(2*r1)<=F(19,20) and F(19,20)**399<=e/8
    and 10/N<=dlo/2 and 400/(dlo*dlo*N)<=e/8,
 'whole_matching_true_source_face_Chebyshev': 240/L<=e and 2*h+1/N<=3*e/4
    and 120*h<=dlo and s/(1-e)<10 and chi-face==F(1,500000)
    and gap>0 and bad<=1-ret,
 'weighted_capture_and_same_tau_payments_T8': nu>0 and nu-1/N>0 and eta>0 and K>=4
    and cap<e and 1/(1-(h/r0)**2)<=1+e and ct*L>=4 and tau<e,
 'positive_roots_full_cost_and_history_sizes': L**2/8>=x and L/4>gamma
    and sig*L/2>=2 and root<slack and 1+r1<6 and rho>20 and N>=2,
 'correct_conservative_S16_tau_quartic': shift==[F(29,6),F(29,2),F(13,2),F(13,12),F(1,16)]
    and all(v>0 for v in shift),
}
paths=[SOURCE/n for n in ('README.md','primary_registration.json','primary.py','primary_result.json',
 'threshold_registration.json','threshold.py','threshold_result.json','manifest.json')]
paths += [PROOF,AUDIT/'full_review.md',AUDIT/'local_uniform.md',AUDIT/'MANIFEST.json']
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
pr=json.loads((SOURCE/'primary_registration.json').read_text());primary=json.loads((SOURCE/'primary_result.json').read_text())
tr=json.loads((SOURCE/'threshold_registration.json').read_text());threshold=json.loads((SOURCE/'threshold_result.json').read_text())
checks['original_frozen_proof_and_receipt_identities']=(hashes[str(PROOF)]==reg['source_proof_sha256']==pr['proof_sha256']
 and hashes[str(SOURCE/'primary_registration.json')]==primary['registration_sha256']
 and hashes[str(SOURCE/'primary.py')]==primary['code_sha256']
 and hashes[str(SOURCE/'primary_result.json')]==tr['primary_result_sha256']
 and hashes[str(SOURCE/'threshold_registration.json')]==threshold['registration_sha256']
 and hashes[str(SOURCE/'threshold.py')]==threshold['script_sha256'])
assert all(checks.values()),checks
values={'N12':N,'rho_lower':rho,'mu_lower':mulo,'mu_upper':muhi,'dmu_lower':dlo,'dmu_upper':dhi,
 'transport_mean_upper':transport_hi,'transport_source_mean_upper':transport_hi/(1-e),
 'baseline_mean_upper':mean,'baseline_failure_upper':4*var/e**2,'count_failure_upper':4/(e*e*N),
 'target_cap_failure_upper':F(19,20)**399,'profile_failure_upper':400/(dlo*dlo*N),
 'source_face_bad_upper':bad,'weighted_exponent_lower':K,'cap_T8_upper':cap,'tau_rate':ct,
 'tau_T8_upper':tau,'full_root_T8_upper':root,'mass_margin':(1-e)*masslo-mass,
 'paid_gamma_margin':(1-e)*Gamma-gamma,'price_margin':price-pricehi,'scalar_slack':slack}
out={'status':'PASS independent threshold conditions; analytic acceptance in review.md','pid':os.getpid(),
 'python':platform.python_version(),'author_primary8_and_threshold56_executed':False,
 'recipe':'exact Fraction; N12; single T8; log32 intervals; one fixed registered tuple and L0',
 'own_registration_sha256':hashlib.sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest(),
 'source_sha256':hashes,'checks':checks,'exact_bounds':{k:str(v) for k,v in values.items()},
 'display_only':{k:float(v) for k,v in values.items()}}
(HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','exact_bounds')},indent=2))
