"""Independent positive-integral Gamma mixture diagnostics; no rectangle claim."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import mpmath as mp
HERE=Path(__file__).resolve().parent

def one(n):
 started=time.monotonic();answers=[]
 for prec in [60,90]:
  with mp.workdps(prec):
   mu=1+mp.lambertw(mp.exp(-1));m=mu*n-1;theta=1+mp.mpf(n)**(-mp.mpf(1)/4);a=m/theta;eps=mp.log(n)/mp.sqrt(n);sd=theta*mp.sqrt(a);sc=[theta,theta*(1-2/mp.sqrt(a)),theta*(1+1/mp.sqrt(a))];lg=mp.loggamma(a)
   def density(z,s):
    t=m+sd*z
    if t<=0:return mp.mpf(0)
    return mp.exp(mp.log(sd)+(a-1)*mp.log(t)-t/s-a*mp.log(s)-lg)
   def cdf(z,index=None):
    weights=[1-eps,eps/3,2*eps/3] if index is None else [mp.mpf(1),mp.mpf(0),mp.mpf(0)]
    fun=lambda zz:mp.fsum(weights[i]*density(zz,sc[i]) for i in range(3) if weights[i])
    return mp.quad(fun,[-20,-10,-5,0,z])
   z=mp.findroot(lambda zz:cdf(zz)-mp.mpf('.5'),(-mp.mpf('.1'),mp.mpf('.1')),tol=mp.mpf(10)**(-prec+15));z0=mp.findroot(lambda zz:cdf(zz,0)-mp.mpf('.5'),(-mp.mpf('.1'),mp.mpf('.1')),tol=mp.mpf(10)**(-prec+15));shift=sd*(z-z0)/mp.log(n)
   answers.append(mp.nstr(shift,42))
 assert answers[0]==answers[1],answers
 return dict(n=n,precision_pair=[60,90],saved_digits=42,median_shift_over_log_n=answers[1],integration_domain_standardized=[-20,'median'],method='positive standardized Gamma density integral; omitted left tail is not part of an interval certificate',elapsed=time.monotonic()-started)
if __name__=='__main__':
 start=time.monotonic();wl=json.loads((HERE/'gamma_obstruction_wolfram.json').read_text());hs=[r['n'] for r in wl['results']]
 with ProcessPoolExecutor(max_workers=12) as p:res=list(p.map(one,hs))
 with mp.workdps(70):
  errors=[abs(mp.mpf(a['median_shift_over_log_n'])-mp.mpf(b['median_shift_over_log_n'].split('`')[0])) for a,b in zip(res,wl['results'])];err=max(errors);assert err<mp.mpf('1e-38');error=mp.nstr(err,12)
 (HERE/'gamma_check.json').write_text(json.dumps(dict(status='complete',workers=12,certified=False,scope='Abstract rank spectrum formula diagnostics only',cross_software_error=error,results=res,elapsed=time.monotonic()-start),indent=2)+'\n');print(json.dumps(dict(status='complete',cases=len(res),cross_software_error=error,elapsed=time.monotonic()-start)))
