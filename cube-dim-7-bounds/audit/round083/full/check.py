#!/usr/bin/env python3
"""Bounded independent R083 sensitivity run; no author imports/execution."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json, math, os, sys
from scipy.integrate import quad
import scipy

HERE=Path(__file__).resolve().parent
reg=json.loads((HERE/'REGISTRATION.json').read_text())
assert all(sha256(Path(p).read_bytes()).hexdigest()==v for p,v in reg['sources'].items())
def log20(y):
    z=(y-1)/(y+1)
    low=2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)),F())
    return low,low+2*z**41/(41*(1-z*z))

e=F(1,100000); s=F(9,2)
qlo,qhi=log20(F(18,11)); muhi=F(7,2)-F(11,2)*qlo
masslo,_=log20(F(4499549,1000001))
_,pricehi=log20(F(1000,387))
oldprice,_=log20(F(5,2))
P0=1-2*e; Ulo=P0*(1-e)/(1+e); Uhi=1/Ulo
M=2*(1+e)*Uhi
b=2*(1-e)*(1-e)*Ulo-2*e
paid=(1-e)*F(1503,1000)*b/M
values={
    'mu_upper':muhi, 'mean_slack':F(7914,10000)-muhi,
    'inverse_mean_slack':F(79147,100000)-F(79144,100000)/(1-e),
    'rate_gap':F(1583,2000)-F(79147,100000),
    'variance_slack':2-F(27,16)/(1-e),
    'target_mass_slack':(1-e)*masslo-F(1503,1000),
    'paid_deficit_slack':paid-F(1501,1000),
    'price_slack':F(19,20)-pricehi,
    'tau_rate':F(613,2000)-F(1583,2000)*F(387,1000),
    'weighted_cap_rate_L':F(22221,100000)/2/F(19998)/(1+e)*F(387,1000),
    'old81_price_slack':oldprice-F(229,250),
    'old81_separator_L64':F(11,20)-(F(48,125)+F(5,64)),
    'old80_separator_L64':F(11,20)-(F(28029,100000)+F(4,64)),
    'fresh_log_derivative':F(3873,88),
    'fresh_J1_derivative':F(9,8000)+F(72,11),
}
assert all(v>0 for v in values.values())
assert values['fresh_log_derivative']<90 and values['fresh_J1_derivative']<40
assert qhi<F(1,2) and values['rate_gap']==F(3,100000)

def formula(w,q,t,h,eps):
    rho=t*q; R=4.5*q-h; a=(1-eps)*q; g=R-a
    m=max(abs(w)-a,0); v=max(g-m,0)
    Q=lambda u:-math.log1p(-u/(2*R))
    FF=lambda u:2*R*Q(u)-u
    return (1+(FF(g)+FF(v)-m)/(2*rho),
            R*(Q(g)**2+Q(v)**2)/(2*rho),
            R*(Q(g)**3+Q(v)**3)/(3*rho))

records=[]; maxerr=0.; maxdiff=[0.,0.,0.]
t=600.; h=.0002; eps=.00001
for q in [.25,.55,.95]:
    rho=t*q; R=4.5*q-h; a=(1-eps)*q; g=R-a
    assert 1<=R<=4.5 and rho>=16 and eps*q<=h/4.5 and 0<=g<=7*R/9
    for tag,w in [('full',0.),('edge',a+.4*g),('central+',R+h),('central-',-R-1.5*h)]:
        assert abs(w)<=R+2*h
        left=max(-rho,w-rho-a); right=min(rho,w+rho+a)
        points=sorted(set([left,right]+[v for v in [-rho+g,0.,rho-g] if left<v<right]))
        direct=[]
        for j in range(3):
            def integrand(z):
                ell=min(1.,(R+rho+a-abs(z))/(2*R))
                phi=-math.log(ell)
                return (phi**j)/ell/(2*rho)
            direct.append(sum(quad(integrand,l,u,epsabs=1e-13,epsrel=1e-13)[0] for l,u in zip(points,points[1:])))
        calc=formula(w,q,t,h,eps)
        err=max(abs(d-c) for d,c in zip(direct,calc)); maxerr=max(maxerr,err)
        dq=1e-6
        plus=formula(w,q+dq,t,h,eps); minus=formula(w,q-dq,t,h,eps)
        diffs=[rho*abs(math.log(plus[0])-math.log(minus[0]))/(2*dq)]+[rho*abs(plus[j]-minus[j])/(2*dq) for j in [1,2]]
        maxdiff=[max(x,y) for x,y in zip(maxdiff,diffs)]
        records.append({'q':q,'regime':tag,'w':w,'max_integral_error':err,'rho_scaled_diagnostics':diffs})
assert maxerr<1e-11
assert all(sha256(Path(p).read_bytes()).hexdigest()==v for p,v in reg['sources'].items())
receipt={'status':'PASS','pid':os.getpid(),'python':sys.executable,'scipy':scipy.__version__,
         'registration_sha256':sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest(),
         'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'rational_values':{k:str(v) for k,v in values.items()},
         'display_only':{k:float(v) for k,v in values.items()},
         'coordinate_cases':records,'max_integral_error':maxerr,'max_rho_scaled_diagnostics':maxdiff,
         'source_hashes_rechecked':True,
         'limits':'Finite quadrature and derivative diagnostics support consistency only. Analytic proof supplies uniformity. No explicit threshold, no author program run, no old baseline replay.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','pid','max_integral_error','max_rho_scaled_diagnostics','source_hashes_rechecked']}))
