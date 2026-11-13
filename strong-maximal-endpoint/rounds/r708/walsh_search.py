"""Exact mixture model for higher sign-parity modes, evaluated in floating point."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import sys,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent
for directory in ['r700','r692','r706','r707']:sys.path.insert(0,str(HERE.parent/directory))
import joint,refine,absolute
from stable_gamma import Special
MU=absolute.MU;P=(MU-1)/MU;T=1-2*P

def run(job):
    n,k=job;start=time.time();L=np.log(n);epsilon=1/np.sqrt(L)
    assert 4*epsilon*(k-1)<27/4 # Full k-action algebraic objective strictly concave.
    _,_,eta=absolute.normal(n);cbase=brentq(lambda c:absolute.derivative(c,eta),.05,4);pbase=sp.expit(-cbase)
    special=Special(64)
    for mod in [joint,refine]:mod.MU=cbase;mod.P=pbase;mod.sp=special
    seed={r['n']:r for r in json.loads((HERE.parent/'r700/refined.json').read_text())}
    ref=refine.run(seed[n]);beta=1+ref['a']*np.sqrt(L/n);w=(1+ref['a']**2)*L/2+np.log(L)+ref['b'];qbase=ref['t1']-w;rate=beta/n
    lo,hi=int(binom.ppf(1e-16,n-k,pbase)),int(binom.isf(1e-16,n-k,pbase))
    j=np.arange(lo,hi+1);wt=binom.pmf(j,n-k,pbase);wt/=sum(wt)
    sides=[]
    for sign in [1,-1]:
        pp=brentq(lambda p:np.log((1-p)/p)-1/(1-p)-2*sign*epsilon*(1-2*p)**(k-1),1e-14,.5,xtol=1e-14)
        cc=np.log((1-pp)/pp);jj=np.arange(k+1);ww=binom.pmf(jj,k,pp)
        shifts=(cbase*j[:,None]+cc*jj[None,:]+sign*epsilon*(-1.)**jj[None,:]).ravel()
        weights=(wt[:,None]*ww[None,:]).ravel()
        def residual(u):
            q=qbase+u;xx=q+w-shifts;ell=-rate*xx-n*np.log1p(-rate)
            return float(weights@(np.expm1(ell)-np.exp(ell)*special.gammainc(n,(1-rate)*np.maximum(xx,0))
                    +special.gammainc(n,np.maximum(xx,0))-special.gammainc(n,np.maximum(q-shifts,0))))
        u=brentq(residual,-2*k-2,2*k+2,xtol=1e-10)
        theory=k*(1+pp*np.log((1-pp)/pp)-MU)+sign*epsilon*(1-2*pp)**k
        sides.append(dict(sign=sign,p=pp,u=u,u_limit=theory,residual=abs(residual(u))))
    gain=sp.logsumexp([r['u'] for r in sides])-np.log(2)-np.log(np.cosh(epsilon))
    kgain=(1-T**(2*k))/2-2*P*(1-P)**2*k*T**(2*k-2)
    return dict(n=n,mode_order=k,epsilon=epsilon,rows=sides,scaled_log_gain=float(L*gain),predicted_scaled_limit=float(-kgain),
        candidate_log_gain=float(gain),max_row_error=max(abs(r['u']-r['u_limit']) for r in sides),
        algebraic_global_concavity_margin=float(27/4-4*epsilon*(k-1)),seconds=time.time()-start)
if __name__=='__main__':
    jobs=[(n,k) for n in [64,1024,16384,262144,1048576] for k in [2,3,4]];out=[]
    with ProcessPoolExecutor(12) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            item=future.result();out.append(item);print(json.dumps(item),flush=True)
            (HERE/'walsh_search.json').write_text(json.dumps(sorted(out,key=lambda r:(r['n'],r['mode_order'])),indent=2)+'\n')
