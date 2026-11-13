"""Optimal one-dimensional scalar profiles at fixed c=mu, using quadrature."""
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import evaluator, MU,P

def run(config):
    n,order,*extras=config
    upper_sd = extras[0] if extras else 14
    mean=n*(1+MU*P);sd=np.sqrt(n*(1+MU*MU*P*(1-P)))
    low,high=int(binom.ppf(1e-15,n,P)),int(binom.isf(1e-15,n,P))
    ks=np.arange(low,high+1);wt=binom.pmf(ks,n,P);wt/=sum(wt)
    logwt=np.log(wt)
    def logdensity(t):
        x=np.maximum(np.atleast_1d(t)[:,None]-MU*ks,1e-300)
        v=x/n-1
        with np.errstate(divide='ignore', invalid='ignore'):
            logv=np.log1p(v)
            stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5)
            logpdf=-.5*np.log(2*np.pi*n)-logv+n*(logv-v)-stir
            logpdf=np.where(np.isfinite(logpdf),logpdf,-np.inf)
        return sp.logsumexp(logpdf+logwt,axis=1)
    def cdf(t):
        return wt@sp.gammainc(n,np.maximum(t-MU*ks,0))
    nodes,weights=leggauss(order)
    upper=mean+upper_sd*sd
    def profile(t0,details=False):
        logf0=logdensity(t0)[0]
        def logr(t):
            return t-t0+logdensity(t)[0]-logf0
        t1=brentq(lambda t:logr(t)-np.log(n),t0,t0+max(100,4*np.log(n)))
        ts=t1+(upper-t1)*(nodes+1)/2
        logf=logdensity(ts)
        target=ts-t0+logf-logf0
        # log Phi'(exp(y)) = (n-2)log(1+y)+log(n+y).
        y=np.maximum((target-np.log(n))/(n-1),0)
        for _ in range(8):
            y-=( (n-2)*np.log1p(y)+np.log(n+y)-target)/((n-2)/(1+y)+1/(n+y))
        pdf=np.exp(logf)
        factor=(upper-t1)/2*weights
        extra=factor@(np.expm1(y)*pdf)
        residual=extra-cdf(t0)
        cost= -np.expm1(-(t1-t0)) + (factor@(np.exp(y)*(1+y)/(n+y)*pdf))/np.exp(logf0)
        score=t0-np.log(cost)
        return (score,t1,residual) if details else residual
    t0=brentq(profile,mean-8*sd,mean,xtol=1e-9)
    score,t1,resid=profile(t0,True)
    baseline,_,_,_=evaluator(n,1e-16)
    base=baseline(0,0)[0]+MU*n
    return dict(n=n,order=order,t0=t0,plateau_width=t1-t0,score=score,
                gain=score-base,average_residual=resid,upper_standard_deviations=upper_sd,
                binomial_terms=len(ks))
