"""Interval evaluation of analytic Chernoff bounds; no distribution tail truncation."""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from mpmath import mp, iv
import numpy as np
from scipy.optimize import brentq

def bound(row, precision=80):
    mp.dps=100
    iv.dps=precision
    root=mp.findroot(lambda c:c-1-mp.exp(-c),1.28)
    low=mp.floor(root*mp.mpf(10)**70)/mp.mpf(10)**70
    high=low+mp.mpf(10)**-70
    left,right=iv.mpf(mp.nstr(low,100)),iv.mpf(mp.nstr(high,100))
    assert left-1-iv.exp(-left)<0 and right-1-iv.exp(-right)>0
    c=iv.mpf([mp.nstr(low,100),mp.nstr(high,100)])
    p=1/(1+iv.exp(c))
    n=row['n'];t0=iv.mpf(str(row['t0']))
    # One positive summand is a rigorous lower bound for the full density.
    j=int(np.floor(n/(1+np.exp(float(root)))))
    x=t0-c*j
    assert x>0
    lg=iv.loggamma
    logf=(lg(n+1)-lg(j+1)-lg(n-j+1)+j*iv.log(p)+(n-j)*iv.log(1-p)
          -x+(n-1)*iv.log(x)-lg(n))
    logD=(-t0-logf+n*iv.log(2*(1-p))-lg(n))/(n-1)
    D=iv.exp(logD)
    T=n*c+row['upper_standard_deviations']*iv.sqrt(n*c)
    cn=float(root);pn=1/(1+np.exp(cn));tn=float(T.mid)
    theta=brentq(lambda s:n/(1-s)+n*cn*pn*np.exp(cn*s)/(1-pn+pn*np.exp(cn*s))-tn,0,0.99)
    th=iv.mpf(str(theta))
    assert 0<D.b<th and th<1
    logmgf=-n*iv.log(1-th)+n*iv.log(1-p+p*iv.exp(c*th))
    mean_bound=logmgf-(th-D.b)*T
    budget_bound=mean_bound-logf.a
    assert mean_bound.b < -50 and budget_bound.b < -50
    return dict(n=n,t0=str(row['t0']),upper_sd=row['upper_standard_deviations'],
                selected_binomial_index=j,theta=str(theta),D_interval=str(D),
                log_mean_tail_interval=str(mean_bound),
                log_scaled_budget_tail_interval=str(budget_bound),
                log10_mean_bound=float((mean_bound.b/iv.log(10)).b),
                log10_scaled_budget_bound=float((budget_bound.b/iv.log(10)).b),
                interval_precision=precision,mu_bracket=[mp.nstr(low,100),mp.nstr(high,100)])
if __name__=='__main__':
    parent=Path(__file__).resolve().parents[1]
    rows=json.loads((parent/'r693/extended.json').read_text())
    with ProcessPoolExecutor(6) as pool:
        out=list(pool.map(bound,rows))
    Path(__file__).with_name('bounds.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in out:print(row['n'],row['log10_mean_bound'],row['log10_scaled_budget_bound'])
