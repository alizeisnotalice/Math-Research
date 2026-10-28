#!/usr/bin/env python3
"""One independent bounded audit check, not the author baseline."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json, math, os, sys
from scipy.integrate import quad
import scipy

HERE=Path(__file__).resolve().parent
reg=json.loads((HERE/'REGISTRATION.json').read_text())
src=Path(reg['source'])
assert sha256(src.read_bytes()).hexdigest()==reg['source_sha256']

def log24(y):
    assert y>1
    z=(y-1)/(y+1)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(24)),F())
    return lo,lo+2*z**49/(49*(1-z*z))

def expminus_half(N):
    return sum((F(-1,2)**j/math.factorial(j) for j in range(N+1)),F())

alo,ahi=expminus_half(25),expminus_half(24)
slo,shi=1/(2*ahi-1),1/(2*alo-1)
clo,_=log24(1/(3-4*alo))
_,chi=log24(1/(3-4*ahi))
log85lo,log85hi=log24(F(8,5))
log53lo,log53hi=log24(F(5,3))
exact={
 'exp_minus_half_low':alo,'exp_minus_half_high':ahi,
 'sstar_low':slo,'sstar_high':shi,'cstar_low':clo,'cstar_high':chi,
 'endpoint4_derivative_numerator_low':1-2*log85hi,
 'endpoint5_derivative_numerator_high':1-2*log53lo,
 'fresh_J_prime_bound':F(81,2),
 'fresh_log_J_prime_bound':F(243,4),
 'fresh_J1_prime_bound':9+F(9,5000),
 'fresh_J1_value':F(9,5),'fresh_J2_value':F(18,25),
 'cap_count_uniform_lower_at_limit':F(1,5)-F(6,5)*F(1,12),
 'true_variance_uniform_upper_at_limit':5/(1-F(1,12)),
}
assert F(3,5)<alo<ahi<F(5,8)
assert F(469348449,100000000)<slo<shi<F(469348450,100000000)
assert F(55533956,100000000)<clo<chi<F(55533957,100000000)
assert exact['endpoint4_derivative_numerator_low']>0
assert exact['endpoint5_derivative_numerator_high']<0
assert exact['fresh_log_J_prime_bound']<120 and exact['fresh_J1_prime_bound']<60
assert exact['cap_count_uniform_lower_at_limit']==F(1,10)
assert exact['true_variance_uniform_upper_at_limit']<10

def profiles(s,w,q,t,h,eps):
    rho=t*q; R=s*q-h; a=(1-eps)*q; g=R-a
    m=max(abs(w)-a,0); v=max(g-m,0)
    Q=lambda u:-math.log1p(-u/(2*R))
    FF=lambda u:2*R*Q(u)-u
    return (1+(FF(g)+FF(v)-m)/(2*rho),
            R*(Q(g)**2+Q(v)**2)/(2*rho),
            R*(Q(g)**3+Q(v)**3)/(3*rho))

records=[]; maxerr=0.; maxdiff=[0.,0.,0.]
t,h,eps=600.,.0002,.00001
for s in [4.,4.5,5.]:
    for q in [.3,.8]:
        rho=t*q; R=s*q-h; a=(1-eps)*q; g=R-a
        assert q>=1/s and 1<=R<=5 and rho>=20 and eps*q<=h/s and 0<=g/R<=1-1/s
        for tag,w in [('full',0.),('edge',a+.6*g),('central+',R+h),('central-',-R-1.5*h)]:
            assert abs(w)<=min(s,R+2*h)
            left=max(-rho,w-rho-a); right=min(rho,w+rho+a)
            points=sorted(set([left,right]+[v for v in [-rho+g,0.,rho-g] if left<v<right]))
            direct=[]
            for j in range(3):
                def fun(z):
                    ell=min(1.,(R+rho+a-abs(z))/(2*R)); phi=-math.log(ell)
                    return phi**j/ell/(2*rho)
                direct.append(sum(quad(fun,l,u,epsabs=1e-13,epsrel=1e-13)[0] for l,u in zip(points,points[1:])))
            f=profiles(s,w,q,t,h,eps)
            err=max(abs(d-v) for d,v in zip(direct,f)); maxerr=max(maxerr,err)
            dq=1e-6; plus=profiles(s,w,q+dq,t,h,eps); minus=profiles(s,w,q-dq,t,h,eps)
            diagnostics=[rho*abs(math.log(plus[0])-math.log(minus[0]))/(2*dq)]+[rho*abs(plus[j]-minus[j])/(2*dq) for j in [1,2]]
            maxdiff=[max(x,y) for x,y in zip(maxdiff,diagnostics)]
            records.append({'s':s,'q':q,'regime':tag,'w':w,'max_integral_error':err,'rho_scaled_derivative_diagnostics':diagnostics})
assert maxerr<1e-11
assert sha256(src.read_bytes()).hexdigest()==reg['source_sha256']
receipt={'status':'PASS','pid':os.getpid(),'python':sys.executable,'scipy':scipy.__version__,
         'registration_sha256':sha256((HERE/'REGISTRATION.json').read_bytes()).hexdigest(),
         'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'source_sha256_rechecked':reg['source_sha256'],'exact_values':{k:str(v) for k,v in exact.items()},
         'display_only':{k:float(v) for k,v in exact.items()},'cases':records,
         'max_integral_error':maxerr,'max_rho_scaled_diagnostics':maxdiff,
         'limits':'Bounded arithmetic and numerical diagnostics only. Uniform theorem and quantifiers require analytic review; no dimension threshold or author baseline replay.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','pid','max_integral_error','max_rho_scaled_diagnostics','source_sha256_rechecked']}))
