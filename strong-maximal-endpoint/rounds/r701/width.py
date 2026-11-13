"""One-dimensional width elimination at fixed critical slope."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import MU,P,evaluator
def run(job):
    n,a=job;L=np.log(n);beta=1+a*np.sqrt(L/n);gamma=1-beta/n
    ev,_,_,_=evaluator(n,1e-16)
    k=np.arange(int(binom.ppf(1e-15,n,P)),int(binom.isf(1e-15,n,P))+1)
    weights=binom.pmf(k,n,P);weights/=sum(weights);lw=np.log(weights)
    cut=n/beta-1;v=(cut-n)/n;stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5)
    logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+n*(v-np.log1p(v))+np.log(sp.gammaincc(n,cut))
    base=(1+a*a)*L/2+np.log(L)
    def ratio(b):
        w=base+b
        logu=-w+np.log(np.expm1(logI))
        logK=np.logaddexp(0,logu)
        q=MU*n+ev(a,b)[0]+logK;t1=q+w
        xx=np.maximum(q-MU*k,1e-300);v=xx/n-1
        with np.errstate(divide='ignore',invalid='ignore'):
            lv=np.log1p(v);ld=-.5*np.log(2*np.pi*n)-lv+n*(lv-v)-stir
            ld=np.where(np.isfinite(ld),ld,-np.inf)
            logf=sp.logsumexp(lw+ld)
            logA=sp.logsumexp(lw-beta*(t1-MU*k)/n-n*np.log(gamma)+np.log(sp.gammaincc(n,gamma*np.maximum(t1-MU*k,0))))
        return logu+np.log(n)+logf-np.log(beta)-logA
    low=max(-8,-base+.01);high=8
    b=brentq(ratio,low,high,xtol=1e-9)
    return dict(n=n,a=float(a),b=float(b),width=float(base+b),
                corrected_constant=float(logI-.5*(1+a*a)*L),
                width_log_ratio_residual=float(ratio(b)),
                log_ratio_derivative_fd=float((ratio(b+.01)-ratio(b-.01))/.02))
if __name__=='__main__':
    jobs=[(n,a) for n in [64,256,4096,65536,1048576,4194304] for a in [.4,.8,np.sqrt(MU),1.6,2.]]
    with ProcessPoolExecutor(12) as pool:out=list(pool.map(run,jobs))
    Path('widths.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in out:
        if abs(r['a']-np.sqrt(MU))<1e-10:print(r,flush=True)
