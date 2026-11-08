import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.special import gammainc,gammaincc,logsumexp
from scipy.stats import binom
from scipy.optimize import brentq
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
def run(j):
 n=j['N'];t=j['shift'];L=math.log(n);rho=math.sqrt(L/n);c=MU+t*rho;p=1/(1+math.exp(c));mean=n*(1+c*p);lo=max(0.,n*MU-7*math.sqrt(n*L));hi=n*MU+7*math.sqrt(n*L)
 radius=math.ceil(math.sqrt(220*n));k=np.arange(max(0,math.floor(n*p)-radius),min(n,math.ceil(n*p)+radius)+1);weights=np.exp(binom.logpmf(k,n,p));low=np.maximum(lo-c*k,0);high=np.maximum(hi-c*k,0);lp=float(weights@gammainc(n,low));up=float(weights@gammaincc(n,high));tail=lp+up;omission=2*math.exp(-440) if len(k)<n+1 else 0.;bound=2*n**-6
 assert c>0 and c<2 and abs(mean-n*MU)<=math.sqrt(n*L) and tail+omission<=bound
 return dict(j,c=c,mean=mean,W_lower=lo,W_upper=hi,W_length=hi-lo,tail=tail,tail_prior_omission_bound=omission,Chernoff_bound=bound,mean_loss_cap2_bound=2*bound,repair_log_cost_bound=-n*math.log1p(-2*bound),mass_residual=float(weights.sum()-1),scope='gamma-binomial floating tail; analytic Chernoff bound independent of floating tail')
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);jobs=[dict(N=n,shift=t) for n in [2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576] for t in [0.,-math.sqrt(MU)/2,-.5,.5]];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'window72.json').write_text(json.dumps(dict(mu=MU,seconds=time.monotonic()-start,workers=12,rows=rows),indent=2)+'\n');print('PASS',len(rows),'wide-window tail checks')
