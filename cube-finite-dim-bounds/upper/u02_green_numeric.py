"""U02: cube killed Green family; A=Delta, source=0, unit domain (-1,1)^n.

Finite deterministic n=3,4 scalar optimizations. Float quadrature errors are
estimates, not interval certificates. Independent mpmath validation is separate.
"""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import erfc

SPLIT = 0.12
IMAGE_M = 4
SPECTRAL_K = 40
ODD = np.arange(1,2*SPECTRAL_K,2,dtype=float)
EIG = (math.pi*ODD/2)**2

def heat_image(t,x):
    if t == 0: return 0.0
    return sum(math.exp(-(x-4*m)**2/(4*t))-
               math.exp(-(x-(4*m+2))**2/(4*t))
               for m in range(-IMAGE_M,IMAGE_M+1))/math.sqrt(4*math.pi*t)

def heat_spectral(t,x):
    return float(np.dot(np.exp(-EIG*t),np.cos(ODD*math.pi*x/2)))

def survival_image(t):
    if t == 0: return 1.0
    return 1-2*sum((-1)**k*erfc((2*k+1)/(2*math.sqrt(t)))
                   for k in range(IMAGE_M+1))

def survival_spectral(t):
    return float(np.dot(np.exp(-EIG*t),4/math.pi*((-1.0)**np.arange(SPECTRAL_K))/ODD))

def integrate_split(fun_lo,fun_hi,n):
    lo,el=quad(lambda t:fun_lo(t)**n,0,SPLIT,epsabs=2e-12,epsrel=2e-12,limit=200)
    hi,eh=quad(lambda t:fun_hi(t)**n,SPLIT,np.inf,epsabs=2e-12,epsrel=2e-12,limit=200)
    return lo+hi,el+eh

def corner(n,a):
    return integrate_split(lambda t:heat_image(t,a),lambda t:heat_spectral(t,a),n)

def main():
    rows=[]
    for n in (3,4):
        mass,mass_error=integrate_split(survival_image,survival_spectral,n)
        opt=minimize_scalar(lambda a:mass/(2**n*a**n*corner(n,a)[0]),
                            bounds=(0.05,0.95),method='bounded',
                            options={'xatol':2e-11,'maxiter':100})
        a=float(opt.x); g,gerr=corner(n,a); r=1/a
        rows.append({'n':n,'a':a,'R':r,'unit_domain_mass':mass,
                     'unit_corner_green':g,'cost':mass/(2**n*a**n*g),
                     'domain_mass_R':r*r*mass,
                     'corner_green_R':r**(2-n)*g,
                     'mass_quad_error_estimate':mass_error,
                     'green_quad_error_estimate':gerr,
                     'optimizer_success':bool(opt.success),'optimizer_nfev':int(opt.nfev)})
    checks=[]
    for t in (0.08,0.12,0.2):
        for x in (0.0,0.4,0.8):
            checks.append({'t':t,'x':x,'heat_absolute_difference':abs(heat_image(t,x)-heat_spectral(t,x)),
                           'survival_absolute_difference':abs(survival_image(t)-survival_spectral(t))})
    result={'generator':'Delta; coordinate variance 2t','rows':rows,
            'expression_checks':checks,'precision_status':'floating point numerical evidence',
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
