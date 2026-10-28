#!/usr/bin/env python3
"""Independent R082 new threshold audit; never executes the author 54 gates.

Weaker N12 and single positive Taylor T8 terms suffice for all sensitive
base margins. S16 tau-tail polynomial identity is checked separately.
This is a single fixed L0, no nearby threshold or parameter search.
"""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import hashlib,json,os,platform

HERE=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/upper/round082/threshold_compute')
PROOF=Path('/Users/zhengzhihao/Desktop/T/cube-r078-theory-bernstein-sparse-20261028/cube-dim-7-bounds/upper/round081/general_s/PROOF.md')
L=F(2**24);x=F(2,5)*L;N=L**12/F(factorial(12));h=1/L**2;e=F(1,10000)
r0=1+F(1,100000);r1=F(99,25)-F(1,100000);q1=F(99,100)
rho=N/(4*x);mu=576*x*x/N+48*x/N;var=2304*x*x/N
P0=1-2*e;Uhi=(1+e)/(P0*(1-e));M=2*(1+e)*Uhi
nu=F(1,4)-F(5,4)*e;eta=F(1,198);K=(nu-1/N)*eta*x/(1+e)
gap=F(29,1000)*x-F(1,2);retbad=2*x/gap**2
T8=lambda y:y**8/F(factorial(8))
cap=M/T8(K);tau=r1*L**2*4096/T8(L/50)
root=1/N+1/T8(F(3,5)*L)+8192/T8(F(3,10)*L)
yc=[-F(5,2),-F(3,2),-F(1,4),F(1,12),F(1,16)]
uc=[sum((yc[j]*comb(j,i)*4**(j-i) for j in range(i,5)),F(0)) for i in range(5)]
checks={
 'explicit_local_and_geometry_domain': N>=300 and x>=2 and x/N<1 and x*N*N>1
    and 2*h<F(1,100000) and h<=F(1,1000) and r0-2*h>1 and 1-h>0
    and F(3,4)-h>0 and 1/(x*N*N)<1 and q1/(x*N*N)<=h/4
    and N/x>r1 and rho>4*r1 and rho>=12 and r1<4 and 4*r1>5,
 'five_event_probability_envelopes': mu<=e/2 and 4*var/e**2<=e/8
    and 4/(e*e*N)<=e/8 and 1+x/(2*N)-1/(2*r1)<=F(9,10)
    and F(9,10)**299<=e/8 and 4/N<=F(1,200) and F(640000)/N<=e/8,
 'whole_matching_and_retention': 80/L<=e and 6/rho<=e and 6/rho<=F(1,2)
    and 2*h+1/N<=3*e/4 and 80*h<=F(1,100)
    and F(67,100)/(1-e)<F(671,1000) and 4/(3*(1-e))<2
    and gap>0 and retbad<=F(1,1000),
 'independent_weighted_cap_and_tau_T8': nu-1/N>0 and K>=4 and cap<e
    and 1/(1-(h/r0)**2)<=1+e and L/50>=4 and tau<e,
 'positive_root_complete_cost_and_finite_history': L**2/8>=x and L/4>F(13,10)
    and F(3,10)*L>=2 and root<F(1,200)
    and N>=2 and rho>5 and 1/N**3-4/(8*N**4)>0,
 'correct_S16_tau_tail_quartic_identity': uc==[F(53,6),F(33,2),F(27,4),F(13,12),F(1,16)]
    and all(v>0 for v in uc),
 'README_cubic_explanation_fails_at_y4': sum(yc[j]*4**j for j in range(4))<0,
}
reg=json.loads((SOURCE/'registration.json').read_text());receipt=json.loads((SOURCE/'result.json').read_text())
sources=[SOURCE/n for n in ('registration.json','README.md','verify.py','result.json')]+[PROOF]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
checks['original_registration_script_and_R081_identity']=(F(reg['L0'])==L
 and hashes[str(PROOF)]==reg['proof_sha256']==receipt['proof_sha256']
 and hashes[str(SOURCE/'registration.json')]==receipt['registration_sha256']
 and hashes[str(SOURCE/'verify.py')]==receipt['script_sha256'])
assert all(checks.values()),checks
values={'L0':L,'N12':N,'x0':x,'h0':h,'rho_lower':rho,'mean_upper':mu,
 'log_failure_upper':4*var/e**2,'count_failure_upper':4/(e*e*N),
 'profile_failure_upper':F(640000)/N,'source_face_bad_upper':retbad,
 'weighted_cap_exponent_lower':K,'cap_T8_upper':cap,'tau_T8_upper':tau,
 'full_root_T8_upper':root}
out={'status':'PASS independent fixed-threshold arithmetic; all-tail/operator proof in review.md',
 'pid':os.getpid(),'python':platform.python_version(),'registration':'single specified L0=2^24, independent audit',
 'method':'exact Fraction; N12; single positive eighth Taylor term; exact polynomial coefficients',
 'author_54_reexecuted':False,'checks':checks,'source_sha256':hashes,
 'bounds':{k:str(v) for k,v in values.items()},'display_only':{k:float(v) for k,v in values.items()},
 'tau_quartic_u_coefficients':[str(v) for v in uc]}
(HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','bounds')},indent=2))
