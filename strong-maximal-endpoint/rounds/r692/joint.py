"""Deterministic joint optimization; finite bounds define the numerical search."""
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq, minimize

MU = float(1+sp.lambertw(1/np.e).real)
P = 1/(1+np.exp(MU))
def evaluator(n, tail):
    L = np.log(float(n))
    low, high = int(binom.ppf(tail,n,P)), int(binom.isf(tail,n,P))
    k = np.arange(low,high+1)
    weights = binom.pmf(k,n,P)
    omitted = float(binom.cdf(low-1,n,P)+binom.sf(high,n,P))
    weights /= weights.sum()
    count = [0]
    def score(a,b,old=False):
        count[0] += 1
        beta = 1+a*np.sqrt(L/n)
        w = (1+a*a)*L/2+np.log(L)+b
        if w <= 0: return -1e20, 1e20
        rate = beta/n
        r = n/beta-1
        v = (r-n)/n
        stir = 1/(12*n)-1/(360*n**3)+1/(1260*n**5)
        logI = (np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+
                n*(v-np.log1p(v))+np.log(sp.gammaincc(n,r)))
        logK = np.logaddexp(np.log(-np.expm1(-w)),logI-w)
        def residual(offset):
            q = MU*n+offset
            x = np.maximum(q+w-MU*k,0)
            ell = -rate*(q+w-MU*k)-n*np.log1p(-rate)
            return weights @ (np.expm1(ell)-np.exp(ell)*sp.gammainc(n,(1-rate)*x)
                    +sp.gammainc(n,x)-sp.gammainc(n,np.maximum(q-MU*k,0)))
        offset = brentq(residual,-MU*n, MU*n+10,xtol=2e-10)
        # Return centered scores to avoid subtracting two large linear terms.
        return offset-logK,float(abs(residual(offset)))
    return score, count, omitted, len(k)

def run(n):
    score,count,omitted,terms = evaluator(n,1e-12)
    base=score(0,0)[0]
    results=[]
    for start in ([0.5,0.5],[1,1.5],[1.5,0]):
        res=minimize(lambda x: -(score(*x)[0]-base),start,
                     method='Nelder-Mead',bounds=[(0,2),(-2,5)],
                     options={'xatol':2e-5,'fatol':2e-8,'maxiter':320})
        results.append(res)
    best=min(results,key=lambda x:x.fun)
    refined,_,tail2,terms2=evaluator(n,1e-16)
    refined_base=refined(0,0)[0]
    gain=refined(*best.x)[0]-refined_base
    previous=refined(np.sqrt(MU),np.log(np.log(n)))[0]-refined_base
    return dict(n=n,a=float(best.x[0]),b=float(best.x[1]),
                gain=gain,previous_gain=previous,tail_stability=abs(gain+best.fun),
                omitted_initial=omitted,omitted_refined=tail2,terms=terms,
                terms_refined=terms2,evaluations=count[0],
                all_starts_success=all(x.success for x in results),
                start_gain_spread=float(max(-x.fun for x in results)-min(-x.fun for x in results)),
                root_residual=refined(*best.x)[1])
if __name__=='__main__':
    ns=[64,128,256,512,1024,4096,16384,65536,262144,1048576,2097152,4194304]
    with ProcessPoolExecutor(12) as pool:
        rows=list(pool.map(run,ns))
    Path(__file__).with_name('results.json').write_text(json.dumps(rows,indent=2)+'\n')
    for row in rows: print(row,flush=True)
