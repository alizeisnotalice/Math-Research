import json
from pathlib import Path
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
from joint import MU,P,evaluator
rows=json.loads(Path('results.json').read_text())
for row in rows:
    n,a,b=row['n'],row['a'],row['b']
    L=np.log(float(n));beta=1+a*np.sqrt(L/n);w=(1+a*a)*L/2+np.log(L)+b
    k=np.arange(int(binom.ppf(1e-16,n,P)),int(binom.isf(1e-16,n,P))+1)
    wt=binom.pmf(k,n,P);wt/=sum(wt)
    rate=beta/n
    def terms(offset):
        q=MU*n+offset
        x=np.maximum(q+w-MU*k,0);ell=-rate*(q+w-MU*k)-n*np.log1p(-rate)
        A=wt@(np.exp(ell)*sp.gammaincc(n,(1-rate)*x))
        residual=wt@(np.expm1(ell)-np.exp(ell)*sp.gammainc(n,(1-rate)*x)
                        +sp.gammainc(n,x)-sp.gammainc(n,np.maximum(q-MU*k,0)))
        xx=q-MU*k
        pos=xx>0
        density=wt[pos]@np.exp((n-1)*np.log(xx[pos])-xx[pos]-sp.gammaln(n))
        return residual,A,density
    root=brentq(lambda x:terms(x)[0],-MU*n,MU*n+10,xtol=2e-10)
    _,A,f=terms(root)
    r=n/beta-1;v=(r-n)/n
    logI=np.log(n/beta)+.5*np.log(2*np.pi/n)+1/(12*n)-1/(360*n**3)+1/(1260*n**5)+n*(v-np.log1p(v))+np.log(sp.gammaincc(n,r))
    u=np.exp(-w)*np.expm1(logI)
    row['width_derivative']=-beta*A/(beta*A+n*f)+u/(1+u)
    row['width_stationarity_ratio']=u*n*f/(beta*A)
    ev,_,_,_=evaluator(n,1e-16);step=1e-3
    row['slope_derivative_fd']=(ev(a+step,b)[0]-ev(a-step,b)[0])/(2*step)
Path('diagnostics.json').write_text(json.dumps(rows,indent=2)+'\n')
print(max(abs(r['width_derivative']) for r in rows))
