"""Two-sign physical-input rows; rectangles may cross both components.
Floating point search, not a finite-n global-optimization certificate.
"""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import sys,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq,minimize_scalar
HERE=Path(__file__).resolve().parent
for dirname in ['r700','r692','r706','r707']:sys.path.insert(0,str(HERE.parent/dirname))
import joint,refine,absolute
from stable_gamma import Special
MU=absolute.MU;PSTAR=(MU-1)/MU;CURV=2*PSTAR**2*(1-PSTAR)

def m(p):return 1 if p==0 else 1+p*np.log((1-p)/p)
def best_limit(delta,sign):
    if -2-2*sign*delta>=0:p=.5
    else:p=brentq(lambda p:np.log((1-p)/p)-1/(1-p)-2*sign*delta,1e-14,.5,xtol=1e-14)
    return p,m(p)-MU+(1-2*p)*sign*delta

def run(job):
    n,alpha=job;start=time.time();L=np.log(n);delta=alpha/np.sqrt(L)
    z,H,eta=absolute.normal(n);cbase=brentq(lambda c:absolute.derivative(c,eta),.05,4);pbase=sp.expit(-cbase)
    special=Special(64)
    for mod in [joint,refine]:mod.MU=cbase;mod.P=pbase;mod.sp=special
    old={r['n']:r for r in json.loads((HERE.parent/'r700/refined.json').read_text())}
    ref=refine.run(old[n]);beta=1+ref['a']*np.sqrt(L/n);w=(1+ref['a']**2)*L/2+np.log(L)+ref['b'];qbase=ref['t1']-w
    lo,hi=int(binom.ppf(1e-16,n-1,pbase)),int(binom.isf(1e-16,n-1,pbase))
    k=np.arange(lo,hi+1);weights=binom.pmf(k,n-1,pbase);weights/=sum(weights)
    omitted=float(binom.cdf(lo-1,n-1,pbase)+binom.sf(hi,n-1,pbase));rate=beta/n
    def receiver(p,sign):
        c=0 if p==0 else np.log((1-p)/p)
        shifts=np.concatenate([cbase*k+sign*delta,cbase*k+c-sign*delta])
        wt=np.concatenate([(1-p)*weights,p*weights])
        def residual(u):
            qq=qbase+u;xx=qq+w-shifts;ell=-rate*xx-n*np.log1p(-rate)
            P1=special.gammainc(n,(1-rate)*np.maximum(xx,0));P2=special.gammainc(n,np.maximum(xx,0));P0=special.gammainc(n,np.maximum(qq-shifts,0))
            return float(wt@(np.expm1(ell)-np.exp(ell)*P1+P2-P0))
        u=brentq(residual,-5-2*delta,5+2*delta,xtol=1e-10)
        return u,abs(residual(u))
    sides=[]
    for sign in [1,-1]:
        plim,ulim=best_limit(delta,sign)
        opt=minimize_scalar(lambda p:-receiver(p,sign)[0],bounds=(1e-8,.5-1e-8),method='bounded',options={'xatol':2e-8,'maxiter':90})
        candidates=[(float(opt.x),receiver(float(opt.x),sign)),(0,receiver(0,sign)),(.5,receiver(.5,sign)),(plim,receiver(plim,sign))]
        bp,best=max(candidates,key=lambda x:x[1][0])
        sides.append(dict(sign=sign,p=bp,p_limit=plim,u=best[0],u_limit=ulim,residual=best[1],optimizer_success=bool(opt.success),evaluations=int(opt.nfev)))
    loggain=sp.logsumexp([r['u'] for r in sides])-np.log(2)-np.log(np.cosh(delta))
    limitgain=sp.logsumexp([r['u_limit'] for r in sides])-np.log(2)-np.log(np.cosh(delta))
    return dict(n=n,alpha=alpha,delta=delta,cbase=cbase,ramp_success=ref['success'],sides=sides,
       log_score_gain=float(loggain),scaled_log_gain=float(L*loggain),
       algebraic_log_gain=float(limitgain),scaled_algebraic_log_gain=float(L*limitgain),
       predicted_scaled_limit=float(-CURV*alpha*alpha),
       row_to_algebra_error=max(abs(r['u']-r['u_limit']) for r in sides),
       omitted_probability=omitted,binomial_terms=len(k),seconds=time.time()-start)
if __name__=='__main__':
    controls='--controls-only' in sys.argv
    jobs=[(n,a) for n in [64,256,1024,4096,16384,65536,262144,1048576] for a in ([0] if controls else [0,.5,1,2])]
    out=json.loads((HERE/'sign_search.json').read_text()) if controls else []
    with ProcessPoolExecutor(12) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();out.append(r);print(json.dumps(r),flush=True)
            (HERE/'sign_search.json').write_text(json.dumps(sorted(out,key=lambda r:(r['n'],r['alpha'])),indent=2)+'\n')
