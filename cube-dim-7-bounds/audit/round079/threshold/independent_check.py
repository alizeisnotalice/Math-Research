#!/usr/bin/env python3
"""Independent R079 fixed-tuple arithmetic. No producer module is executed.

This directory is independent temporary staging, NOT a Git clone.
Uses N32, 24-term atanh logarithms and P48, unlike producer N64/log64/P64.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib
import json
import os
import platform

ROOT = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
HERE = Path(__file__).resolve().parent

def log_interval(v):
    assert v > 1
    z = (v-1)/(v+1)
    total = sum((2*z**(2*j+1)/Q(2*j+1) for j in range(24)), Q(0))
    return total, total+2*z**49/(49*(1-z*z))

def exp_prefix(y):
    assert y > 0
    return sum((y**j/Q(factorial(j)) for j in range(49)), Q(0))

# Explicit independent fixed inputs. Registration reconciliation is checked below.
L=Q(4096); s=Q(3,2); q0=Q(2,3); q1=Q(89,100)
sigma=Q(11,100); chi=Q(51,1000); D=Q(1536)
delta=Q(1,2**19); eB=eN=eT=delta
eJ=eR=eS=eD=Q(1,2**20); epsA=Q(1,2**16)
ret=Q(999999,1000000); mass=Q(2889,10000)
gamma=Q(7166,25000); price=Q(5827,50000); c=Q(17,100); slack=Q(1,10000)
x=(1-sigma)*L; h=1/(D*L); N=L**32/Q(factorial(32))
r0=1+delta; r1=s*q1-delta
P0=1-eB-6/D; ulo=P0*(1-eN)/(1+eJ); uhi=(1+eN)/(P0*(1-eJ))
M=2*(1+eR)*uhi; b=2*ret*(1-eR)*ulo-eS-eT; G=mass*b/M; budget=(1-eD)*G
qm,qp=log_interval(Q(6,5)); masslo,_=log_interval(r1/r0)
mu1hi=Q(1,2)-Q(5,2)*qm
mu2lo=1-5*qp-Q(5,2)*qp*qp; mu2hi=1-5*qm-Q(5,2)*qm*qm
meanlo=(1-1/N)*(mu2lo-h*(mu2hi+Q(1,25)))
meanhi=mu2hi+1/(250*N*q0)
m=Q(2211,50000)/(1-eJ); p=Q(183,1000); v=Q(531,100000)
gap=(chi-m)*x-p; H=gap*gap/(2*(v*x+p*gap/3))
eta=(1-q1)/(2*q1); cS=(eta-chi)*(1-sigma); cT=sigma/2-chi*(1-sigma)
_,source_loghi=log_interval(q1/q0); _,pricehi=log_interval(1/(1-sigma))
cap=Q(1025,1024)*source_loghi*r1*D*L/exp_prefix(cS*L)
# sqrt(L0)=64 exactly. The proof uses the L^(3/2) envelope after this base.
tau=r1*D*L*64/exp_prefix(cT*L)
scalar=1/N+1/exp_prefix(sigma*L)+128/exp_prefix(sigma*L/2)
ML=Q(3,4)*L**2/(q0*q0*N)+Q(1,2)*L/(q0*N)
VL=Q(3,2)*L**2/(q0*q0*N)

groups={
 'full_active_geometry': (s*q0==1 and q0<q1<1 and 1<r0<r1<Q(3,2)
    and 2*h<delta and h<=Q(1,2**20) and q1/(x*N*N)<=h/s
    and 1/(x*N*N)<Q(1,2) and N*q0/x>4*r1 and N/x>r1
    and x/N<1 and x*N*N>1 and 2<=x<=L and N>=9000 and x/N<=Q(1,1000)),
 'target_six_event_budgets': (ML<=eB/2 and 4*VL/eB**2<=epsA/8
    and 4/(eN**2*N)<=epsA/8 and 1+x/(2*N)-1/(2*r1)<=Q(3,4)
    and N>=75 and Q(3,4)**74<=epsA/8 and mu1hi<Q(221,5000)
    and Q(324000000)/N<=epsA/8 and 0<mu2lo<mu2hi
    and meanlo>Q(528,100000) and meanhi<Q(5295,1000000)
    and Q(5760000)/N<=epsA/8),
 'whole_matching_profiles': (2*h+1/N<=3*eN/4 and x/(3*q0*N)<=eJ<=Q(1,1024)
    and h<=Q(1,100000) and h/5<Q(1,1000000)
    and Q(5301,1000000)/(1-eJ)<v and qp<p),
 'uniform_Bernstein': (chi>m and gap>0 and exp_prefix(H)>1000000),
 'radial_mass_and_degrees': (1/(1-(h/r0)**2)<=1+eR and 0<P0<1
    and M>2 and b>0 and (1-epsA)*masslo>mass and 0<gamma<budget<G<1),
 'unweighted_cap_and_tau': (cS==Q(961,100000) and cT==cS
    and eta*x/N<=Q(1,1025) and cap<eS and tau<eT
    and cS*L>=4 and cT*L>=4),
 'positive_threshold_and_full_root': ((L/2)**2/2>=L and L/4>gamma
    and pricehi<price and sigma*L/2>=4 and scalar<slack
    and gamma-price-c==slack),
}

reg=json.loads((ROOT/'upper/round079/compute/certificate_registration.json').read_text())
amd=json.loads((ROOT/'upper/round079/compute/certificate_amendment.json').read_text())
strong=json.loads((ROOT/'upper/round079/compute/stronger_registration.json').read_text())
chosen=reg['chosen_parameters']
expected={'L0':L,'s':s,'qmin':q0,'qmax':q1,'sigma':sigma,'chi':chi,'hden':D,
    'delta':delta,'eB':eB,'eN':eN,'eT':eT,'eJ':eJ,'eRad':eR,'eS':eS,'eD':eD,'epsA':epsA}
groups['registration_reconciliation']=(all(Q(chosen[k])==val for k,val in expected.items())
    and Q(amd['target_mass_lower'])==mass and Q(amd['Ret'])==ret
    and Q(strong['gamma'])==gamma and Q(strong['price'])==price
    and Q(strong['c'])==c and Q(strong['slack'])==slack)

paths=['upper/round073/cantelli_frontier/PROOF.md','upper/round078/bernstein_theory/PROOF.md',
 'upper/round079/compute/certificate.py','upper/round079/compute/certificate_registration.json',
 'upper/round079/compute/certificate_amendment.json','upper/round079/compute/certificate_result.json',
 'upper/round079/compute/stronger_registration.json','upper/round079/compute/stronger_guard.py',
 'upper/round079/compute/stronger_guard.json','audit/round073/cantelli_frontier/review.md']
hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
receipt=json.loads((ROOT/'upper/round079/compute/stronger_guard.json').read_text())
groups['preserved_source_receipt_hashes']=(receipt['registration_sha256']==hashes[paths[6]]
 and receipt['certificate_sha256']==hashes[paths[5]] and receipt['script_sha256']==hashes[paths[7]]
 and amd['parent_registration_sha256']==hashes[paths[3]])

assert all(groups.values()), [k for k,vv in groups.items() if not vv]
values={'N32':N,'x0':x,'h0':h,'P0':P0,'Ulo':ulo,'Uhi':uhi,'M':M,'b':b,
 'Gamma':G,'budget':budget,'budget_minus_gamma':budget-gamma,
 'target_mass_lower_enclosure':(1-epsA)*masslo,'mass_margin':(1-epsA)*masslo-mass,
 'mu2_lower':mu2lo,'mu2_upper':mu2hi,'mean2_lower':meanlo,'mean2_upper':meanhi,
 'Bernstein_gap':gap,'Bernstein_H0':H,'P48_H0':exp_prefix(H),'Ret_lower':1-1/exp_prefix(H),
 'eta':eta,'cap_rate':cS,'tau_rate':cT,'cap_upper':cap,'tau_upper':tau,
 'price_upper':pricehi,'price_margin':price-pricehi,'root_remainder_upper':scalar,'slack':slack}
out={'status':'PASS independent arithmetic; analytic acceptance is in review.md',
 'pid':os.getpid(),'python':platform.python_version(),'staging':'independent directory; NOT a clone',
 'method':'Fraction; N32; 24 positive atanh terms plus exact tail; positive P48',
 'producer_programs_executed':False,'groups':groups,'source_sha256':hashes,
 'bounds':{k:str(val) for k,val in values.items()},'display_only':{k:float(val) for k,val in values.items()}}
(HERE/'independent_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'groups':groups,
 'display_only':out['display_only']},indent=2))
