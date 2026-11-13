"""Budget-preserving compact perturbations of the optimized finite ramp."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import MU,P
def run(job):
    row,mom,order=job;n=row['n'];L=np.log(n);a,b=row['a'],row['b']
    beta=1+a*np.sqrt(L/n);w=(1+a*a)*L/2+np.log(L)+b;t1=mom['t1']
    z=n/beta-1;v=(z-n)/n;stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5)
    logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+n*(v-np.log1p(v))+np.log(sp.gammaincc(n,z))
    lower=max(-5,(t1-MU*n)/np.sqrt(n)+.5);upper=5
    nodes,weights=leggauss(order)
    edges=sorted(set([lower,lower+.5,4.,upper]))
    xx=np.concatenate([lo+(hi-lo)*(nodes+1)/2 for lo,hi in zip(edges[:-1],edges[1:])])
    factor=np.concatenate([(hi-lo)/2*weights*np.sqrt(n) for lo,hi in zip(edges[:-1],edges[1:])])
    ts=MU*n+np.sqrt(n)*xx
    chi=np.sin(np.pi/2*np.clip((xx-lower)/.5,0,1))**2*np.cos(np.pi/2*np.clip(xx-4,0,1))**2
    mf,vf=mom['source_center'],mom['source_variance']
    mg,vg=mom['budget_center'],mom['budget_variance']
    poly=.5*np.log(vg/vf)+(xx-mg)**2/(2*vg)-(xx-mf)**2/(2*vf)
    u=ts-t1;y=beta*u/n
    logg=-u+beta*u/n+(n-1)*np.log1p(beta*u/n)-logI
    def cost(kappa):
        dy=chi*(poly-kappa)/n
        ratio=dy+(n-1)*np.log1p(dy/(1+y))
        return factor@(np.exp(logg)*np.expm1(ratio))
    kappa=brentq(cost,float(poly.min()-1),float(poly.max()+1),xtol=1e-13)
    k=np.arange(int(binom.ppf(1e-15,n,P)),int(binom.isf(1e-15,n,P))+1)
    wt=binom.pmf(k,n,P);wt/=sum(wt)
    x=np.maximum(ts[:,None]-MU*k,1e-300);v=x/n-1
    with np.errstate(divide='ignore',invalid='ignore'):
        lv=np.log1p(v);lp=-.5*np.log(2*np.pi*n)-lv+n*(lv-v)-stir
        lp=np.where(np.isfinite(lp),lp,-np.inf)
    logf=sp.logsumexp(lp+np.log(wt),axis=1)
    gain=factor@(np.exp(logf+y)*np.expm1(chi*(poly-kappa)/n))
    return dict(n=n,order=order,kappa=float(kappa),n_mean_gain=float(n*gain),
                normalized_soft_budget_residual=float(cost(kappa)),lower=lower,upper=upper)
if __name__=='__main__':
    rows=json.loads(Path('../r692/results.json').read_text())
    moments={r['n']:r for r in json.loads(Path('../r697/moments.json').read_text())}
    jobs=[(r,moments[r['n']],o) for r in rows if r['n'] in [256,1024,16384,262144,1048576,4194304] for o in [256,512]]
    with ProcessPoolExecutor(6) as pool:out=list(pool.map(run,jobs))
    Path('compact.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in out:print(r,flush=True)
