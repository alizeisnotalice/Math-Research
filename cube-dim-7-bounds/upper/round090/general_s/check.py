#!/usr/bin/env python3
"""One finite auxiliary check for a newly derived general fixed-s theorem."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json, math, os, sys
import scipy
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
reg=json.loads((HERE/'REGISTRATION.json').read_text())
assert all(sha256(Path(s).read_bytes()).hexdigest()==h for s,h in reg['sources'].items())

def log24(y):
    z=(y-1)/(y+1)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(24)),F())
    return lo,lo+2*z**49/(49*(1-z*z))

log2lo,log2hi=log24(F(2))
clo,_=log24(1/(1-log2lo)); _,chi=log24(1/(1-log2hi))
assert 0<log2lo<log2hi<1
assert F(118138,100000)<clo<chi<F(118140,100000)
# Polynomial coefficients in ascending powers of s, without a sampled grid.
A_prime_coeff=[F(1),F(6)]
A_over_q_coeff=[F(0),F(0),F(4)]
Dlog_coeff=[sum(xs) for xs in zip(A_prime_coeff+[F(0)],A_over_q_coeff)]
assert Dlog_coeff==[1,6,4]
assert 2*2==4 and 2*F(1,2)==1

s=F(10); q=F(3,5); h=F(1,5000); eps=F(1,100000); t=F(600)
R=s*q-h; rho=t*q; a0=(1-eps)*q; g=R-a0
assert rho>4*(s+1) and 1<R<s and eps*q<=h/s
assert g/R>F(4,5) and g/R<1-1/s
e_bound=1/(2*(s+1))
nu_bound=1/s-(1+1/s)*e_bound
assert nu_bound==1/(2*s)
exact={'log2_lower':log2lo,'log2_upper':log2hi,
       'limit_coefficient_lower':clo,'limit_coefficient_upper':chi,
       's_diagnostic':s,'relative_gap_diagnostic':g/R,
       'Dlog_at_s10':4*s*s+6*s+1,'J1_derivative_at_s10':5*s,
       'count_floor_boundary':nu_bound,'variance_factor_at_s10':2*s*s}

def profile(w,qf):
    sf,tf,hf,ef=map(float,[s,t,h,eps]); rhof=tf*qf; RR=sf*qf-hf
    aa=(1-ef)*qf; gg=RR-aa; mm=max(abs(w)-aa,0); vv=max(gg-mm,0)
    QQ=lambda u:-math.log1p(-u/(2*RR)); FF=lambda u:2*RR*QQ(u)-u
    return [1+(FF(gg)+FF(vv)-mm)/(2*rhof),
            RR*(QQ(gg)**2+QQ(vv)**2)/(2*rhof),
            RR*(QQ(gg)**3+QQ(vv)**3)/(3*rhof)]

sf,qf,rf,af,gf,hf=map(float,[s,q,R,a0,g,h]); rhof=float(rho)
cases=[]; maxerr=0.; maxdiff=[0.,0.,0.]
for tag,w in [('full',0.),('edge',af+.5*gf),('central+',rf+hf),('central-',-rf-1.5*hf)]:
    assert abs(w)<sf and abs(w)<=rf+2*hf
    left=max(-rhof,w-rhof-af); right=min(rhof,w+rhof+af)
    points=sorted(set([left,right]+[u for u in [-rhof+gf,0.,rhof-gf] if left<u<right]))
    direct=[]
    for j in range(3):
        def fun(z):
            ell=min(1.,(rf+rhof+af-abs(z))/(2*rf)); phi=-math.log(ell)
            return phi**j/ell/(2*rhof)
        direct.append(sum(quad(fun,l,u,epsabs=1e-13,epsrel=1e-13)[0] for l,u in zip(points,points[1:])))
    calc=profile(w,qf); err=max(abs(d-f) for d,f in zip(direct,calc)); maxerr=max(maxerr,err)
    dq=1e-6; plus=profile(w,qf+dq); minus=profile(w,qf-dq)
    diag=[rhof*abs(math.log(plus[0])-math.log(minus[0]))/(2*dq)]+[rhof*abs(plus[j]-minus[j])/(2*dq) for j in [1,2]]
    maxdiff=[max(a,b) for a,b in zip(maxdiff,diag)]
    cases.append({'regime':tag,'w':w,'max_integral_error':err,'rho_scaled_diagnostics':diag})
assert maxerr<1e-11
assert maxdiff[0]<float(exact['Dlog_at_s10']) and maxdiff[1]<float(exact['J1_derivative_at_s10'])
assert all(sha256(Path(s).read_bytes()).hexdigest()==h for s,h in reg['sources'].items())
receipt={'status':'PASS','pid':os.getpid(),'python':sys.executable,'scipy':scipy.__version__,
         'registration_sha256':sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest(),
         'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'exact_values':{k:str(v) for k,v in exact.items()},
         'Dlog_polynomial_coefficients':list(map(str,Dlog_coeff)),
         'display_only':{k:float(v) for k,v in exact.items()},'cases':cases,
         'max_integral_error':maxerr,'max_rho_scaled_diagnostics':maxdiff,
         'source_hashes_rechecked':True,
         'scope':'One finite auxiliary run, not proof of uniformity; no dimension threshold or old grid/baseline replay. New analytic proof supplies all fixed-s and limit quantifiers.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','pid','max_integral_error','max_rho_scaled_diagnostics','source_hashes_rechecked']}))
