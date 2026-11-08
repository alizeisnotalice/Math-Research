import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.special import logsumexp
from scipy.optimize import brentq
from scipy.stats import binom
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
def posterior(n,c,s):
 last=min(n,math.ceil(s/c)-1)
 def lr(k):
  g=s-c*k
  return -math.inf if g<=c else math.log(n-k)-math.log(k+1)+(n-1)*math.log1p(-c/g)
 lo,hi=0,last
 while lo<hi:
  k=(lo+hi)//2
  if lr(k)>0:lo=k+1
  else:hi=k
 m=lo;radius=math.ceil(math.sqrt(110*(n+2)));k=np.arange(max(0,m-radius),min(last,m+radius)+1);p=(m+.5)/(n+1);logw=binom.logpmf(k,n,p)-(k-m)*math.log(p/(1-p))+(n-1)*np.log1p(-c*(k-m)/(s-c*m));w=np.exp(logw-logsumexp(logw));return k,w

def run(j):
 n=j['N'];L=math.log(n);rho=math.sqrt(L/n);D=1.5 if j['D_type']=='constant' else L**(-.25);ref=MU-math.sqrt(MU)/2*rho;act=ref+D*rho;z=brentq(lambda x:.5*math.erfc(x/math.sqrt(2))-math.sqrt(MU/n)*x,0,math.sqrt(L)+2);s=n*MU if j['regime']=='mean' else n*MU-math.sqrt(MU*n)*z;k,w=posterior(n,act,s);target=1/(n*L*L);tails=np.cumsum(w[::-1])[::-1];idx=int(np.flatnonzero(tails<=target)[0]);cut=int(k[idx]);ka,wa=posterior(n,ref,s);at=float(tails[idx]);rt=float(wa[ka>=cut].sum());p=1/(1+math.exp(MU))
 return dict(j,D=D,c_reference=ref,c_actual=act,S=s,cut=cut,actual_removed=at,target_removed=target,reference_removed=rt,reference_loss_over_L_over_N=rt*n/L,action_scalar_log_penalty_model=.5*p*D*D*L,posterior_truncation_bound=(n+1)*math.exp(-220),scope='conditional count diagnostic at one S; not a physical full-profile counterexample or certificate')
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=2**k,D_type=kind,regime=reg) for k in [6,8,10,12,14,16,18,20] for kind in ['constant','vanishing'] for reg in ['mean','critical']];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'posterior73.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,workers=12,rows=rows),indent=2)+'\n');print('PASS',len(rows),'conditional posterior probes')
