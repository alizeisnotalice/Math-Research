"""Critical quantile prediction and joint-stationarity moment diagnostics."""
import json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'r692'))
from joint import MU,P
def run(job):
    row,mom=job;n=row['n'];L=np.log(n);beta=mom['beta'];t1=mom['t1'];rate=1-beta/n
    z=brentq(lambda x:sp.log_ndtr(-x)-.5*np.log(MU)-np.log(x)+.5*L,.01,10)
    pred=np.sqrt(MU)*z/np.sqrt(L)
    cut=n/beta-1;v=(cut-n)/n;stir=1/(12*n)-1/(360*n**3)+1/(1260*n**5)
    logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+stir+n*(v-np.log1p(v))+np.log(sp.gammaincc(n,cut))
    tilt=beta/n;pt=P*np.exp(tilt*MU)/(1-P+P*np.exp(tilt*MU))
    k=np.arange(int(binom.ppf(1e-15,n,pt)),int(binom.isf(1e-15,n,pt))+1)
    wt=binom.pmf(k,n,pt);wt/=sum(wt)
    x=np.maximum(rate*(t1-MU*k),0);Q=sp.gammaincc(n,x)
    with np.errstate(divide='ignore'):
        h=np.exp(n*np.log(x)-x-sp.gammaln(n))
    mx=(wt@((MU*k+n/rate-t1)*Q+h/rate))/(wt@Q)
    egx=(mom['budget_center']*np.sqrt(n)+MU*n)-t1
    rhs=(egx-1)/(-np.expm1(-logI))
    return dict(n=n,a=row['a'],a_quantile=pred,a_difference=row['a']-pred,
                scaled_a_difference=(row['a']-pred)*np.sqrt(L),
                tilted_source_x=mx,budget_x=egx,
                beta_stationarity_residual=mx-rhs,
                budget_minus_tilted_source_mean=egx-mx,
                predicted_mean_difference=1-mx*np.exp(-logI),
                logI=logI)
if __name__=='__main__':
    rows=json.loads(Path('../r692/results.json').read_text())
    moments={r['n']:r for r in json.loads(Path('../r697/moments.json').read_text())}
    with ProcessPoolExecutor(12) as pool:out=list(pool.map(run,[(r,moments[r['n']]) for r in rows]))
    Path('stationary.json').write_text(json.dumps(out,indent=2)+'\n')
    for x in out:print(x,flush=True)
