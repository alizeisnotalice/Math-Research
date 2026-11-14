#!/usr/bin/env python3
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import json,time,sys,platform,resource,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import scipy
from scipy import special,integrate
HERE=Path(__file__).resolve().parent
NS=[16,64,256,1024,4096,16384]
MU=1+float(special.lambertw(np.exp(-1)).real);Q=1/(1+np.exp(MU))
def calc(task):
 n,tol=task;j=np.arange(n+1,dtype=float);m=MU*n;b=n**.25;t=n**(-.75)
 lw=special.gammaln(n+1)-special.gammaln(j+1)-special.gammaln(n-j+1)+j*np.log(Q)+(n-j)*np.log1p(-Q)
 norm=float(special.logsumexp(lw));w=np.exp(lw-norm);lg=special.gammaln(n)
 def loggam(x):
  pos=x>0;out=np.full_like(x,-np.inf);out[pos]=(n-1)*np.log(x[pos])-x[pos]-lg;return out
 def pdf(x):return float(np.exp(special.logsumexp(lw-norm+loggam(x-MU*j))))
 den=special.logsumexp(lw-norm+loggam(m-MU*j))
 num=special.logsumexp(lw-norm-(n-j)*t+loggam(m-MU*j-(n-j)*t))
 sc=float(np.exp(num-den));su=float(np.exp(n*np.log1p((1-Q)*np.expm1(-t))))
 lo=np.maximum(m-2*b-MU*j,0);hi=np.maximum(m+2*b-MU*j,0)
 cp=float(np.dot(w,special.gammainc(n,hi)-special.gammainc(n,lo)))
 cq,err=integrate.quad(pdf,m-2*b,m+2*b,epsabs=tol,epsrel=tol,limit=100)
 return dict(n=n,quad_tolerance=tol,conditional_survival=sc,unconditional_survival=su,central_probability_cdf=cp,central_probability_positive_integral=cq,quad_error_estimate=err,cdf_integral_difference=cp-cq,tv_lower_expression=1-sc-cp-su,tv_lower_nonnegative=max(0,1-sc-cp-su),binomial_log_normalization_residual=norm,central_density=float(np.exp(den)),worker_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
if __name__=='__main__':
 start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool: rows=list(pool.map(calc,[(n,tol) for n in NS for tol in [1e-9,1e-11]]))
 out=dict(status='complete',scope='Floating evaluation of an analytic TV lower bound; not interval certification or C_n counterexample',dimensions=NS,workers=12,blas_threads=1,mu=MU,q=Q,versions=dict(python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()),wall_seconds=time.monotonic()-start,results=rows)
 out['pure_exponential_exact_formula_float']=[dict(n=n,conditional_event_probability=-math.expm1((n-1)*math.log1p(-math.log(n)/n)),tv_lower_expression=-math.expm1((n-1)*math.log1p(-math.log(n)/n))-math.log(n)*math.exp((n-2)*math.log(n-2)-(n-2)-math.lgamma(n-1))) for n in [16,64,256,1024,10000,1000000]]
 (HERE/'rank_tv_bounds.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'cases':len(rows),'wall_seconds':out['wall_seconds'],'max_cdf_integral_difference':max(abs(r['cdf_integral_difference']) for r in rows),'values':[[r['n'],r['tv_lower_nonnegative']] for r in rows if r['quad_tolerance']==1e-11]}))
