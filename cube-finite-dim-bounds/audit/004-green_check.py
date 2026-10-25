#!/usr/bin/env python3
"""Fixed a=1/2, n=3 Green cross-check; floating quadrature, not an interval proof."""
import hashlib,json,math
from pathlib import Path
from scipy.integrate import quad
from scipy.special import erfc
N=3;A=0.5;R=2.;SWITCH=0.15;TERMS=12

def heat(t):
    if t==0:return 0.
    if t<=SWITCH:
        return sum(math.exp(-(A+4*k)**2/(4*t))-math.exp(-(A+2+4*k)**2/(4*t)) for k in range(-TERMS,TERMS+1))/math.sqrt(4*math.pi*t)
    return sum(math.exp(-((2*k+1)*math.pi/2)**2*t)*math.cos((2*k+1)*math.pi*A/2) for k in range(TERMS))
def survive(t):
    if t==0:return 1.
    if t<=SWITCH:
        return 1-2*sum((-1)**k*erfc((2*k+1)/(2*math.sqrt(t))) for k in range(TERMS))
    return 4/math.pi*sum((-1)**k/(2*k+1)*math.exp(-((2*k+1)*math.pi/2)**2*t) for k in range(TERMS))

def integrate(fun):
    intervals=[(0,SWITCH),(SWITCH,1),(1,8),(8,math.inf)]
    runs=[quad(fun,lo,hi,epsabs=2e-13,epsrel=2e-13,limit=100) for lo,hi in intervals]
    return sum(v for v,e in runs),sum(e for v,e in runs)
g,ge=integrate(lambda t:heat(t)**N)
m,me=integrate(lambda t:survive(t)**N)
cost=m/(2**N*A**N*g)
result={'scope':'n=3 R=2 a=1/2 only; no parameter optimization','generator':'Delta; variance2t Brownian motion','switch_time_unit_domain':SWITCH,'terms_each_series':TERMS,'unit_domain_green_corner':g,'unit_domain_green_quadrature_error_estimate':ge,'unit_domain_mean_exit':m,'unit_domain_exit_quadrature_error_estimate':me,'radius2_green_corner':R**(2-N)*g,'radius2_mean_exit':R**2*m,'kernel_cost':cost,'status':'numerical quadrature only; errors exclude truncation/roundoff certification','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('004-green_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
