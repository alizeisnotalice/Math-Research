"""Refine approximate optima with analytic width and slope stationary equations."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import root
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import MU,P,evaluator
def run(row):
    n=row['n'];L=np.log(n);ev,_,_,_=evaluator(n,1e-16)
    k=np.arange(int(binom.ppf(1e-15,n,P)),int(binom.isf(1e-15,n,P))+1)
    wt=binom.pmf(k,n,P);wt/=sum(wt)
    stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5)
    def state(par):
        a,b=par;beta=1+a*np.sqrt(L/n);w=(1+a*a)*L/2+np.log(L)+b
        rate=1-beta/n;cut=n/beta-1;v=(cut-n)/n
        logQ=np.log(sp.gammaincc(n,cut))
        logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+n*(v-np.log1p(v))+logQ
        I=np.exp(logI);u=np.exp(-w)*np.expm1(logI);K=1+u
        q=MU*n+ev(a,b)[0]+np.log(K);t1=q+w
        xx=np.maximum(q-MU*k,1e-300);v=xx/n-1
        logf=sp.logsumexp(-.5*np.log(2*np.pi*n)-np.log1p(v)+n*(np.log1p(v)-v)-stir+np.log(wt))
        pt=P*np.exp(beta*MU/n)/(1-P+P*np.exp(beta*MU/n))
        kt=np.arange(int(binom.ppf(1e-15,n,pt)),int(binom.isf(1e-15,n,pt))+1)
        wtt=binom.pmf(kt,n,pt);wtt/=sum(wtt)
        x=np.maximum(rate*(t1-MU*kt),0);Q=sp.gammaincc(n,x)
        with np.errstate(divide='ignore'):
            h=np.exp(n*np.log(x)-x-sp.gammaln(n))
        mass=wtt@Q
        A=np.exp(-beta*t1/n-n*np.log(rate)+n*np.log1p(P*np.expm1(beta*MU/n)))*mass
        mx=(wtt@((MU*kt+n/rate-t1)*Q+h/rate))/mass
        v=(cut-n)/n;tau=np.exp(.5*np.log(n/(2*np.pi))+n*(np.log1p(v)-v)-stir-logQ)
        egx=n*((beta-1)+beta/n)/(rate*beta)+tau/rate
        residual=[(mx-I/(I-1)*(egx-1))/np.sqrt(n),u*n*np.exp(logf)/(beta*A)-1]
        return residual,mx,egx,logI,q,t1
    sol=root(lambda x:state(x)[0],[row['a'],row['b']],tol=1e-9)
    res,mx,egx,logI,q,t1=state(sol.x)
    return dict(n=n,a=float(sol.x[0]),b=float(sol.x[1]),success=bool(sol.success),
                solver_message=str(sol.message),evaluations=int(sol.nfev),
                residual=list(map(float,res)),budget_minus_tilted_mean=float(egx-mx),
                predicted_mean_difference=float(1-mx*np.exp(-logI)),
                gain_change=float(ev(*sol.x)[0]-ev(row['a'],row['b'])[0]),t1=float(t1))
if __name__=='__main__':
    rows=json.loads(Path('../r692/results.json').read_text())
    with ProcessPoolExecutor(12) as pool:out=list(pool.map(run,rows))
    Path('refined.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in out:print(r,flush=True)
