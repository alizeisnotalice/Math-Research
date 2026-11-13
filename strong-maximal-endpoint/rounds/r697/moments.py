"""Exact incomplete-Gamma moment formulas for the soft budget density."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import MU,P,evaluator
def run(row):
    n=row['n'];a,b=row['a'],row['b'];L=np.log(n)
    beta=1+a*np.sqrt(L/n);w=(1+a*a)*L/2+np.log(L)+b
    rate=1-beta/n;z=n/beta-1
    stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5);v=(z-n)/n
    logQ=np.log(sp.gammaincc(n,z))
    logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+n*(v-np.log1p(v))+logQ
    logK=np.logaddexp(np.log(-np.expm1(-w)),logI-w)
    ev,_,_,_=evaluator(n,1e-16);q=MU*n+ev(a,b)[0]+logK;t1=q+w
    # Stable boundary mass z^n exp(-z)/Gamma(n).
    logtau=.5*np.log(n/(2*np.pi))+n*(np.log1p(v)-v)-stir-logQ
    tau=np.exp(logtau)
    m0=t1-n/beta+n/rate
    mg=m0+tau/rate
    vg=(n+tau*(z+1-n)-tau*tau)/rate**2
    k=np.arange(int(binom.ppf(1e-15,n,P)),int(binom.isf(1e-15,n,P))+1)
    wt=binom.pmf(k,n,P);wt/=sum(wt)
    x=np.maximum(t1-MU*k,0);tail=sp.gammaincc(n,x)
    with np.errstate(divide='ignore'):
        h=np.exp(n*np.log(x)-x-sp.gammaln(n))
    off=MU*(k-n*P);mass=wt@tail
    m1=wt@(off*tail+h)/mass
    m2=wt@((off*off+n)*tail+h*(2*off+x+1-n))/mass
    vf=(m2-m1*m1)/n
    mgc=(mg-MU*n)/np.sqrt(n);mfc=m1/np.sqrt(n);vgc=vg/n
    gaussian=.5*(np.log(vgc/vf)+(vf+(mfc-mgc)**2)/vgc-1)
    return dict(n=n,beta=beta,t1=t1,budget_center=mgc,budget_variance=vgc,
                unconditioned_budget_center=(m0-MU*n)/np.sqrt(n),
                removed_budget_probability=float(sp.gammainc(n,z)),
                source_center=mfc,source_variance=vf,moment_gaussian_kl=gaussian)
if __name__=='__main__':
    rows=json.loads((Path(__file__).resolve().parents[1]/'r692/results.json').read_text())
    with ProcessPoolExecutor(12) as pool: out=list(pool.map(run,rows))
    Path(__file__).with_name('moments.json').write_text(json.dumps(out,indent=2)+'\n')
    for x in out:print(x,flush=True)
