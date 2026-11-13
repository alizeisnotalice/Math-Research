"""Run with BLAS threads set to one; independent dimensions use a process pool."""
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy import special as sp
from scipy.optimize import brentq
from scipy.stats import binom

mu = 1 + sp.lambertw(1/np.e).real
prob = 1/(1+np.exp(mu))

def dimension(n):
    L = np.log(float(n))
    lo, hi = binom.ppf([1e-14, 1-1e-14], n, prob).astype(int)
    k = np.arange(lo, hi+1)
    weights = binom.pmf(k,n,prob)
    weights /= weights.sum()
    def score(a, old=False):
        beta = 1+a*np.sqrt(L/n)
        w = (1+a*a)*L/2+(1 if old else 2)*np.log(L)
        rate = beta/n
        r = n/beta-1
        # Stirling expansion cancels large nearly equal log-gamma terms.
        v = (r-n)/n
        stir = (sp.gammaln(n)-(n-.5)*np.log(n)+n-.5*np.log(2*np.pi)
                if n<100 else 1/(12*n)-1/(360*n**3)+1/(1260*n**5))
        logI = (np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+
                n*(v-np.log1p(v))+np.log(sp.gammaincc(n,r)))
        logK = np.logaddexp(np.log(-np.expm1(-w)),logI-w)
        def residual(q):
            x = np.maximum(q+w-mu*k,0)
            ell = -rate*(q+w-mu*k)-n*np.log1p(-rate)
            return weights @ (np.expm1(ell)-np.exp(ell)*sp.gammainc(n,(1-rate)*x)
                              +sp.gammainc(n,x)-sp.gammainc(n,np.maximum(q-mu*k,0)))
        q = brentq(residual,0.,2*mu*n+10,xtol=1e-11)
        return q-logK,abs(residual(q))
    base,_ = score(0,True)
    rows=[]
    for a in np.unique(np.r_[np.linspace(0,2,33),np.sqrt(mu)]):
        val,err=score(a)
        rows.append((n,a,val-base,np.sqrt(mu)*a-a*a/2,err))
    return rows

if __name__=='__main__':
    ns=sorted(set([2**j for j in range(2,21)]+[10**j for j in range(2,7)]))
    with ProcessPoolExecutor(12) as pool:
        groups=list(pool.map(dimension,ns))
    with Path(__file__).with_name('scores.csv').open('w') as out:
        writer=csv.writer(out);writer.writerow(['n','a','gain','limit','root_residual'])
        for group in groups: writer.writerows(group)
    for group in groups:
        print(group[0][0], max(group,key=lambda z:z[2])[1:3],
              min(group,key=lambda z:abs(z[1]-np.sqrt(mu)))[2],flush=True)
