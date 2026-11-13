"""Exact grouped Bernoulli mixtures up to explicitly recorded tail truncation."""
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
MU=float(1+sp.lambertw(1/np.e).real)
def run(job):
    n,pattern=job;L=np.log(n);rho=np.sqrt(L/n);c0=MU-np.sqrt(MU)/2*rho
    if pattern=='uniform': groups=[(n,c0)]
    elif pattern=='balanced':groups=[(n//2,c0+.75*rho),(n-n//2,c0-.75*rho)]
    elif pattern=='zero_log':groups=[(int(L),0.),(n-int(L),c0)]
    else:groups=[(int(L),5.),(n-int(L),c0)]
    shifts=np.array([0.]);weights=np.array([1.]);mean=float(n);variance=float(n);omitted=0.
    for count,c in groups:
        p=1/(1+np.exp(c));mean+=count*c*p;variance+=count*c*c*p*(1-p)
        if c==0:continue
        lo=int(binom.ppf(1e-14,count,p));hi=int(binom.isf(1e-14,count,p))
        ks=np.arange(lo,hi+1);ww=binom.pmf(ks,count,p)
        omitted+=float(binom.cdf(lo-1,count,p)+binom.sf(hi,count,p))
        shifts=(shifts[:,None]+c*ks).ravel()
        weights=(weights[:,None]*ww).ravel()
    weights/=weights.sum();lw=np.log(weights)
    def residual(q):
        x=np.maximum(q-shifts,0)
        ell=-(q-shifts)/n-n*np.log1p(-1/n)
        return weights@(np.expm1(ell)-np.exp(ell)*sp.gammainc(n,(1-1/n)*x))
    q=brentq(residual,0,2*mean,xtol=1e-9)
    def density(t):
        x=np.maximum(t-shifts,1e-300);v=x/n-1
        with np.errstate(divide='ignore',invalid='ignore'):
            lv=np.log1p(v);lp=-.5*np.log(2*np.pi*n)-lv+n*(lv-v)-(1/(12*n)-1/(360*n**3)+1/(1260*n**5))
            lp=np.where(np.isfinite(lp),lp,-np.inf)
        return float(np.exp(sp.logsumexp(lw+lp)))
    f=density(q);F=float(weights@sp.gammainc(n,np.maximum(q-shifts,0)));d=mean-q
    return dict(n=n,pattern=pattern,groups=groups,mixture_states=len(weights),
                omitted_probability_bound=omitted,deficit_over_log=(n*MU-mean)/L,
                variance_over_n=variance/n,soft_root=q,root_residual=abs(float(residual(q))),
                n_density_over_log=n*f/L,tail_mass_over_d_n=F/(d/n),
                mills_ratio=f/F/(d/variance),neighbor_density_ratio=density(q+L)/f)
if __name__=='__main__':
    jobs=[(n,p) for n in [64,256,1024,4096,16384] for p in ['uniform','balanced','zero_log','far5_log']]
    with ProcessPoolExecutor(12) as pool:rows=list(pool.map(run,jobs))
    Path('triangular.json').write_text(json.dumps(rows,indent=2)+'\n')
    for row in rows:
        if row['n']==16384:print(row,flush=True)
