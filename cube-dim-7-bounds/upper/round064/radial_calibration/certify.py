#!/usr/bin/env python3
"""Exact fixed-tuple radial budget. No grid, quadrature or MATLAB."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform

HERE = Path(__file__).resolve().parent
LOCAL = Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md')
EXPECTED_LOCAL = 'c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1'
assert hashlib.sha256(LOCAL.read_bytes()).hexdigest() == EXPECTED_LOCAL

s = F(3,2)
qmin, qmax = F(2,3), F(909,1000)
inset, hc = F(1,2**30), F(1,2**64)
rlo, rhi = 1+inset, s*qmax-inset
e = F(1,1024)
sigma, chi = F(91,1000), F(1,20)
P0 = 1-e-F(1,2**50)
Ulo, Uhi = P0*(1-e)/(1+e), (1+e)/(P0*(1-e))
M = 2*(1+e)*Uhi
b = 2*(1-e)*(1-e)*Ulo-2*e
mL = F(1239,4000)
Gamma = mL*b/M
gamma_budget = (1-e)*Gamma
gamma = F(61,200)

# One new target log, bounded by a positive series with a rational tail.
z = (rhi-rlo)/(rhi+rlo)
terms = 64
loglower = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
logupper = loglower+2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
assert (1-e)*loglower > mL

# The odd alternating Taylor polynomial is an upper bound on log(6/5).
u = F(1,5)
Q5 = sum(((-1)**(j+1)*u**j/F(j) for j in range(1,6)),F(0))
assert Q5 == F(34187,187500) < F(1,5)
mu = s*Q5*Q5/(1-e)
var = F(2,3)*s*Q5**3/(1-e)
assert mu < F(999,20000) < chi
assert var < F(1,150)
assert chi-F(999,20000) == F(1,20000)

eta = (1-qmax)/(2*qmax)
assert eta-chi == F(1,18180) > 0
assert sigma/2-chi*(1-sigma) == F(1,20000) > 0
assert F(0)<qmin<qmax<1 and 1<rlo<rhi<F(3,2)
assert rhi == F(2727,2000)-inset
assert 0<P0<1 and M>2 and b>0
assert 0<gamma<gamma_budget<Gamma<1
assert 512*hc == F(1,2**55) < F(1,2**50)
craw = gamma-sigma/(1-sigma)
assert craw == F(37249,181800)
assert craw-F(1,5) == F(889,181800) > 0

fixed = dict(s=s,qmin=qmin,qmax=qmax,inset=inset,h_coefficient=hc,
             sigma=sigma,chi=chi,epsilon_J=e,epsilon_N=e,
             epsilon_Ret=e,epsilon_S=e,epsilon_T=e,epsilon_D=e,
             epsilon_A=e,radial_width_tolerance=e,baseline_log_tolerance=e)
derived = dict(rlo=rlo,rhi=rhi,P0=P0,Ulo=Ulo,Uhi=Uhi,M=M,b=b,
               target_mass_lower=mL,Gamma_lower=Gamma,
               gamma_budget_lower=gamma_budget,gamma_chosen=gamma,
               Q5=Q5,meanrate_upper=mu,meanrate_simple=F(999,20000),
               varrate_upper=var,varrate_simple=F(1,150),
               caprate=eta,source_cap_exponent_gap=eta-chi,
               tau_exponent_gap=sigma/2-chi*(1-sigma),
               source_logmass_upper=(qmax-qmin)/qmin,
               product_bridge_upper=512*hc,c_raw=craw,
               margin_above_one_fifth=craw-F(1,5))
result = {
    'status': 'Exact fixed-tuple certificate for the proposed R064 mechanism; full theorem remains pending independent review.',
    'python': platform.python_version(), 'arithmetic': 'fractions.Fraction',
    'target_log_method': '64 positive atanh terms and an exact geometric tail upper bound',
    'pin_log_method': 'five-term odd alternating Taylor upper bound at 1/5',
    'parameters': {k:str(v) for k,v in fixed.items()},
    'derived_rationals': {k:str(v) for k,v in derived.items()},
    'checks': ['target mass lower includes epsilon_A', 'new relative pin-log mean and variance',
               'positive source-cap and true tau exponents', 'radial degree budget M>2',
               'chosen deficit below complete energy budget', 'complete scalar margin above one fifth'],
    'source_sha256': {'relative_coordinate_lemma':EXPECTED_LOCAL},
    'all_assertions_passed': True,
}
(HERE/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({
    'all_assertions_passed':True, 'gamma_chosen':str(gamma),
    'gamma_budget_display':float(gamma_budget), 'c_raw':str(craw),
    'margin_above_one_fifth':str(craw-F(1,5)),
    'source_cap_exponent_gap':str(eta-chi),
    'tau_exponent_gap':str(sigma/2-chi*(1-sigma)),
    'relative_coordinate_lemma':EXPECTED_LOCAL,
},indent=2))
