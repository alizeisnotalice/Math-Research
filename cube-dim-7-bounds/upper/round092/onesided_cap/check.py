#!/usr/bin/env python3
"""One bounded R092 conditional critical-scaling diagnostic, not an operator test."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
import json, math, os, sys
from scipy.integrate import quad
from scipy.special import exp1
import scipy

HERE=Path(__file__).resolve().parent
reg=json.loads((HERE/'REGISTRATION.json').read_text())
assert all(sha256(Path(s).read_bytes()).hexdigest()==h for s,h in reg['sources'].items())
def log24(y):
    z=(y-1)/(y+1)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(24)),F())
    return lo,lo+2*z**49/(49*(1-z*z))

llo,lhi=log24(F(2)); lam=F(3); kap=F(19,20); alpha=F(1,40)
rates={'lambda_critical_margin_lower':1-lam*(1-llo),
       'source_threshold_gap_lower':kap-lam*(1-llo),
       'tau_over_k_rate':1-kap-alpha}
assert all(v>0 for v in rates.values())
# Exact endpoint checks include a negative w far outside |w|<a0.
rho,r,a0,q=F(100),F(8),F(4,5),F(4,5)
for w in [-r,a0-F(1,1000000)]:
    assert w<a0 and -rho>w-rho-a0
    assert -rho+(1-q)<w+rho+a0

x=2.; Q=.75; n=1000000; t=n/x
integrand=lambda q: 0. if q<=0 else math.exp(-x*(1/q-1)/2)/q
I,errI=quad(integrand,0.,1.,epsabs=2e-13,epsrel=2e-13)
J,errJ=quad(integrand,0.,Q,epsabs=2e-13,epsrel=2e-13)
I_E1=math.exp(x/2)*float(exp1(x/2))
J_E1=math.exp(x/2)*float(exp1(x/(2*Q)))
assert abs(I-I_E1)<1e-12 and abs(J-J_E1)<1e-12
def finite_integrand(q):
    if q<1/t: return 0.
    return math.exp((n-1)*math.log1p(-x*(1/q-1)/(2*n)))/q
In,errn=quad(finite_integrand,1/t,1.,epsabs=2e-13,epsrel=2e-13)
naive=math.log(float(lam))+math.log(Q)+I/2
capcost=J/2
paid=naive-capcost
penalty=math.log(Q)+(I-J)/2
assert penalty<0 and paid<math.log(float(lam))
assert abs(In-I)<1e-5
assert all(sha256(Path(s).read_bytes()).hexdigest()==h for s,h in reg['sources'].items())
receipt={'status':'PASS','pid':os.getpid(),'python':sys.executable,'scipy':scipy.__version__,
         'registration_sha256':sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest(),
         'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'exact_rates':{k:str(v) for k,v in rates.items()},
         'fixed_case':{'x':x,'lambda':str(lam),'qmax':Q,'kappa':str(kap),'alpha':str(alpha),'n_diagnostic':n},
         'B_integral_limit':I,'B_integral_E1':I_E1,'finite_n_B_integral':In,
         'finite_n_minus_limit':In-I,'source_cap_integral':J,'source_cap_E1':J_E1,
         'naive_model_coefficient':naive,'integrated_positive_pin_cap_cost':capcost,
         'paid_conditional_model_coefficient':paid,'log_lambda':math.log(float(lam)),
         'qmax_net_penalty':penalty,
         'quadrature_error_estimates':[errI,errJ,errn],
         'source_hashes_rechecked':True,
         'scope':'One fixed conditional scaling model and one-sided endpoint sanity checks only. Not a Cn bound, no n-dependent-s full-chain theorem or grid.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','pid','B_integral_limit','finite_n_minus_limit','naive_model_coefficient','integrated_positive_pin_cap_cost','paid_conditional_model_coefficient','log_lambda']}))
