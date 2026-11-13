import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from mpmath import mp
def check(row):
    mp.dps=60;n=row['n'];t=mp.mpf(str(row['t0']))
    c=1+mp.lambertw(1/mp.e);p=1/(1+mp.exp(c))
    cf=float(c);pf=float(p)
    width=16*np.sqrt(n*pf*(1-pf))
    lo=max(0,int(np.floor(n*pf-width)));hi=min(n,int(np.ceil(n*pf+width)))
    k=np.arange(lo,hi+1);x=float(t)-cf*k
    weights=binom.pmf(k,n,pf)
    direct=(n-1)*np.log(x)-x-sp.gammaln(n)
    v=x/n-1;lv=np.log1p(v)
    stable=-.5*np.log(2*np.pi*n)-lv+n*(lv-v)-(1/(12*n)-1/(360*n**3)+1/(1260*n**5))
    old=sp.logsumexp(direct+np.log(weights))
    new=sp.logsumexp(stable+np.log(weights))
    weight=mp.binomial(n,lo)*p**lo*(1-p)**(n-lo)
    total=mp.mpf(0);lg=mp.loggamma(n)
    for j in range(lo,hi+1):
        xx=t-c*j
        total+=weight*mp.exp((n-1)*mp.log(xx)-xx-lg)
        weight*=mp.mpf(n-j)/(j+1)*p/(1-p)
    exactlog=mp.log(total)
    return dict(n=n,precision=60,terms=hi-lo+1,reference_log_density=mp.nstr(exactlog,50),
                direct_log_error=float(mp.mpf(old)-exactlog),
                stable_log_error=float(mp.mpf(new)-exactlog))
if __name__=='__main__':
    rows=[x for x in json.loads(Path('rounds/r695/stable.json').read_text()) if x['order']==512]
    with ProcessPoolExecutor(3) as pool:r=list(pool.map(check,rows))
    Path('rounds/r695/density_checks.json').write_text(json.dumps(r,indent=2)+'\n')
    print(r)
