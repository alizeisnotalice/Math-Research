"""Finite-dimensional entropy of the optimized soft budget distribution."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import MU,P,evaluator
def run(row):
    n=row['n'];L=np.log(n);a,b=row['a'],row['b']
    beta=1+a*np.sqrt(L/n);w=(1+a*a)*L/2+np.log(L)+b
    r=n/beta-1;v=(r-n)/n
    stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5)
    logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+n*(v-np.log1p(v))+np.log(sp.gammaincc(n,r))
    logK=np.logaddexp(np.log(-np.expm1(-w)),logI-w)
    score,_,_,_=evaluator(n,1e-16);q=score(a,b)[0]+logK+MU*n
    t1=q+w;upper=MU*n+18*np.sqrt(MU*n)
    k=np.arange(int(binom.ppf(1e-15,n,P)),int(binom.isf(1e-15,n,P))+1)
    wt=binom.pmf(k,n,P);wt/=sum(wt)
    mass=wt@sp.gammaincc(n,np.maximum(t1-MU*k,0))
    out=[]
    for order in [256,512]:
        z,weights=leggauss(order);ts=t1+(upper-t1)*(z+1)/2
        x=np.maximum(ts[:,None]-MU*k,1e-300);v=x/n-1
        with np.errstate(divide='ignore',invalid='ignore'):
            lv=np.log1p(v);lp=-.5*np.log(2*np.pi*n)-lv+n*(lv-v)-stir
            lp=np.where(np.isfinite(lp),lp,-np.inf)
        logf=sp.logsumexp(lp+np.log(wt),axis=1)
        xx=ts-t1
        logg=-xx+beta*xx/n+(n-1)*np.log1p(beta*xx/n)-logI
        kl=(upper-t1)/2*weights@(np.exp(logf)/mass*(logf-np.log(mass)-logg))
        out.append(dict(n=n,order=order,conditional_source_mass=float(mass),kl=float(kl)))
    return out
if __name__=='__main__':
    rows=json.loads(Path('../r692/results.json').read_text())
    with ProcessPoolExecutor(12) as pool:groups=list(pool.map(run,rows))
    out=[x for g in groups for x in g]
    Path('entropy.json').write_text(json.dumps(out,indent=2)+'\n')
    for x in out:
        if x['order']==512:print(x,flush=True)
