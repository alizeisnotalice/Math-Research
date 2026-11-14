"""Central Gamma-binomial quantiles: diagnostics, not uniform certificates."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json,re
import numpy as np
from scipy.special import gammainc,lambertw
from scipy.stats import binom
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent
mu=1+float(lambertw(np.exp(-1)).real)
def run(n):
 j=np.arange(n+1);w=binom.pmf(j,n,(mu-1)/mu);target=.5+1/n
 def cdf(q):return float(np.dot(w,gammainc(n+1,np.maximum(q-mu*j,0))))
 q,info=brentq(lambda x:cdf(x)-target,mu*n-10,mu*n+20,xtol=1e-11,rtol=1e-14,full_output=True)
 return dict(n=n,q=q,offset=q-mu*n,cdf_residual=cdf(q)-target,converged=info.converged,iterations=info.iterations,weight_sum=float(w.sum()),offset_limit=(3+mu)/6)
def wlfloat(s):
 return float(re.sub(r'`[^* ]*','',s).replace('*^','e'))
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,[2**k for k in range(2,17)]))
 wl=json.loads((HERE/'quantile_wolfram.json').read_text());checks=[]
 for n in [2**k for k in range(2,10)]:
  pair=[r for r in wl['results'] if r['n']==n];py=next(r for r in rows if r['n']==n)
  checks.append(dict(n=n,python_wolfram_abs_difference=abs(py['q']-wlfloat(pair[1]['q'])),all_wolfram_converged=all(r['status']=='converged' for r in pair)))
 out=dict(workers=12,rows=rows,wolfram_checks=checks,scope='Finite central quantiles only; not general H_n or C_n certificate',seed='deterministic',scipy_method='binomial pmf and central regularized incomplete gamma; Brent root')
 (HERE/'quantile_python.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(cases=len(rows),max_residual=max(abs(r['cdf_residual']) for r in rows),max_cross_difference=max(r['python_wolfram_abs_difference'] for r in checks)),indent=2))
