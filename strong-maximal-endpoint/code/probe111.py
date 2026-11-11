"""Central Gamma–Binomial quantiles for the explicit equal-action cross witness."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import sys,json,math
import numpy as np
from scipy.special import gammainc,ndtr
from scipy.stats import binom
from scipy.optimize import brentq

def one(spec):
 n,u=spec;mu=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=5e-15);sigma=math.sqrt(mu);c=mu+sigma*u/(2*math.sqrt(n));prob=1/(1+math.exp(c));target=ndtr(u);ks=np.arange(int(binom.ppf(1e-15,n,prob)),int(binom.isf(1e-15,n,prob))+1);ps=binom.pmf(ks,n,prob);tail=binom.cdf(ks[0]-1,n,prob)+binom.sf(ks[-1],n,prob)
 def eq(q):return float(np.dot(ps,gammainc(n,np.maximum(0,q-c*ks))))-target
 q=brentq(eq,0,3*n+10*math.sqrt(n),xtol=1e-10);offset=q-mu*n-sigma*math.sqrt(n)*u;pred=(3-mu)/6*(u*u-1)+(mu-1)/8*u*u;assert abs(eq(q))<1e-9 and tail<1e-13
 return dict(N=n,u=u,c=c,rank=q,offset=offset,predicted_offset=pred,offset_difference=offset-pred,cdf_residual=abs(eq(q)),truncated_probability=tail)
if __name__=='__main__':
 specs=[(2**j,u/4) for j in range(5,19) for u in range(-10,11)]
 with ProcessPoolExecutor(12) as pool:rows=list(pool.map(one,specs))
 Path(sys.argv[1]).write_text(json.dumps(dict(cases=len(rows),scope='Fixed action numerical witnesses at central probabilities',max_CDF_residual=max(x['cdf_residual'] for x in rows),max_omitted_binomial_mass=max(x['truncated_probability'] for x in rows),selected=[x for x in rows if x['N'] in (32,1024,16384,262144) and x['u'] in (-2,-1,0,1,2)]),indent=2,default=float)+'\n')
