"""Positive Fenchel-gap decomposition for lambda=n times the ramp soft budget."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import MU,P,evaluator
def run(job):
    row,order=job;n=row['n'];L=np.log(n);a,b=row['a'],row['b']
    beta=1+a*np.sqrt(L/n);w=(1+a*a)*L/2+np.log(L)+b
    cut=n/beta-1;v=(cut-n)/n;stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5)
    logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+n*(v-np.log1p(v))+np.log(sp.gammaincc(n,cut))
    logu=-w+np.log(np.expm1(logI));logK=np.logaddexp(0,logu)
    ev,_,_,_=evaluator(n,1e-16);q=MU*n+ev(a,b)[0]+logK;t1=q+w
    loglam=np.log(n)-q-w+logI
    k=np.arange(int(binom.ppf(1e-15,n,P)),int(binom.isf(1e-15,n,P))+1)
    wt=binom.pmf(k,n,P);wt/=sum(wt);lw=np.log(wt)
    def logf(t):
        x=np.maximum(np.atleast_1d(t)[:,None]-MU*k,1e-300);v=x/n-1
        with np.errstate(divide='ignore',invalid='ignore'):
            lv=np.log1p(v);lp=-.5*np.log(2*np.pi*n)-lv+n*(lv-v)-stir
            lp=np.where(np.isfinite(lp),lp,-np.inf)
        return sp.logsumexp(lp+lw,axis=1)
    target=lambda t:loglam+t+logf(t)[0]
    t0=brentq(target,q-2,q+2)
    ts=brentq(lambda t:target(t)-np.log(n),q-2,t1+10)
    upper=MU*n+18*np.sqrt(MU*n)
    edges=sorted(set([min(q,t0),q,t0,t1,ts,upper]))
    nodes,weights=leggauss(order);pieces=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        tt=lo+(hi-lo)*(nodes+1)/2;fac=(hi-lo)/2*weights
        lf=logf(tt);ls=loglam+tt+lf
        yr=beta*np.maximum(tt-t1,0)/n
        ys=np.zeros_like(ls);active=ls>np.log(n)
        z=(ls[active]-np.log(n))/(n-1)
        for _ in range(9):
            z-=((n-2)*np.log1p(z)+np.log(n+z)-ls[active])/((n-2)/(1+z)+1/(n+z))
        ys[active]=z
        g=np.zeros_like(tt)
        zero=tt<q;plat=(tt>=q)&(tt<t1);soft=tt>=t1
        mask=zero&(ls>0)&~active;g[mask]=-np.expm1(-ls[mask])
        mask=zero&active;g[mask]=np.exp(ys[mask])*(1-(1+ys[mask])/(n+ys[mask]))
        mask=plat&(ls<0);g[mask]=np.expm1(-ls[mask])
        mask=plat&active;g[mask]=np.exp(-ls[mask])+np.expm1(ys[mask])-np.exp(ys[mask])*(1+ys[mask])/(n+ys[mask])
        mask=soft&active
        cr=(1+yr[mask])/(n+yr[mask]);cs=(1+ys[mask])/(n+ys[mask])
        logr=(n-2)*np.log1p(yr[mask])+np.log(n+yr[mask])-ls[mask]
        g[mask]=np.exp(yr[mask])*(np.expm1(logr)*cr+(n-1)*(yr[mask]-ys[mask])/((n+yr[mask])*(n+ys[mask]))+np.expm1(ys[mask]-yr[mask])*(1-cs))
        assert np.all(active[soft])
        val=float(fac@(n*np.exp(lf)*g))
        label='soft' if (lo+hi)/2>=t1 else 'plateau_or_zero'
        pieces.append((label,val))
    total=sum(v for _,v in pieces);small=sum(v for typ,v in pieces if typ!='soft')
    fraction=np.exp(-w+logI-logK)
    return dict(n=n,order=order,normalized_gap=total,normalized_nonsoft_gap=small,
                scaled_score_bound=float(-L*np.log1p(-fraction*total)),
                soft_budget_fraction=float(fraction),t0=t0,t1=t1,tstar=ts)

if __name__=='__main__':
    rows=json.loads(Path('../r700/refined.json').read_text())
    rows=[r for r in rows if r['n'] in [64,1024,16384,262144,1048576,4194304]]
    with ProcessPoolExecutor(6) as pool:out=list(pool.map(run,[(r,o) for r in rows for o in [256,512]]))
    Path('gap.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in out:
        if r['order']==512: print(r,flush=True)
