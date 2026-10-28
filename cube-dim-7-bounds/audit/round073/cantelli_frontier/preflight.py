"""Independent R073 sensitivity check, not a frozen full-theorem audit."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json, os

L=Q(8192); sigma=Q(25,256); chi=Q(63,1280)
x=(1-sigma)*L; qmin=Q(2,3); qmax=Q(11569,12800)
h=1/(768*L); delta=Q(1,2**19); rlo=1+delta; rhi=Q(3,2)*qmax-delta
e=Q(1,2**20); en=Q(1,2**19); ea=Q(1,65536)
N=L**32/factorial(32)
mean=Q(2211,50000)/(1-e); pin=Q(183,1000); var=Q(531,100000)
g=(chi-mean)*x-pin
ret=g*g/(var*x+g*g)
p=1-en-Q(1,128)
ulo=p*(1-en)/(1+e); uhi=(1+en)/(p*(1-e))
Gamma=Q(3043,10000)*(2*Q(9717,10000)*(1-e)*ulo-e-en)/(2*(1+e)*uhi)

def log8(v):
    z=(v-1)/(v+1)
    lo=sum(2*z**(2*j+1)/(2*j+1) for j in range(8))
    return lo,lo+2*z**17/(17*(1-z*z))
def p40(y):
    return sum(y**j/factorial(j) for j in range(41))
eta=(1-qmax)/(2*qmax)
cS=(eta-chi)*(1-sigma); cT=sigma/2-chi*(1-sigma)
cap=Q(1025,1024)*log8(qmax/qmin)[1]*rhi*768*L/p40(cS*L)
tau=rhi*768*L*91/p40(cT*L)
price=sum(sigma**j/j for j in range(1,9))+sigma**9/(9*(1-sigma))
checks={
 'qmax_and_x_are_distinct':qmax!=1-sigma and cS!=cT,
 'wide_h_trim_and_count':delta/h==12 and 2*h+1/N<3*en/4,
 'fixed_retention_and_monotone_domain':ret>Q(9717,10000) and chi>mean and x>pin/(chi-mean),
 'tight_mass_and_deficit':(1-ea)*log8(rhi/rlo)[0]>Q(3043,10000) and (1-e)*Gamma>Q(291,1000),
 'small_finite_n_cap_correction':eta*x/N<Q(1,1025),
 'accurate_cap_and_separate_tau':cap<e and tau<en and min(cS*L,cT*L)>4,
 'full_scalar_price':price<Q(257,2500) and Q(291,1000)-Q(257,2500)-Q(47,250)==Q(1,5000),
}
assert all(checks.values()),checks
out={'status':'PASS','pid':os.getpid(),'checks':checks,'Ret_base_lower':str(ret),
     'Gamma':str(Gamma),'cap_bound':str(cap),'tau_bound':str(tau),
     'scope':'Unfrozen tuple preflight only. N32, 8-term positive atanh log with tail, P40; no author helper or old program rerun.'}
Path(__file__).with_name('preflight.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','pid':out['pid'],'groups':len(checks)}))
