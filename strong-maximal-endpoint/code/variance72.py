import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln,logsumexp,gammainc

def run(j):
 n=j['N'];x=j['Ac'];fac=j['refinement'];ks=np.arange(n+1);logb=gammaln(n+1)-gammaln(ks+1)-gammaln(n-ks+1)
 def angular(r):
  if r<=x:return 1.
  u=r/x;k=ks[ks<u];g=u-k;lw=logb[:len(k)]+(n-1)*np.log(g);p=np.exp(lw-logsumexp(lw));return float(p@((g/u)**2+(n+1)/(n-1)*k*(n-k)/u**2))
 def fn(r):return .5*math.exp(-r)*r*r*angular(r)
 cuts=[0.]+[x*k for k in range(1,n+1) if x*k<80.]+[80.]
 vals=[quad(fn,a,b,epsabs=2e-11*fac/max(1,len(cuts)-1),epsrel=2e-10*fac,limit=100) for a,b in zip(cuts,cuts[1:])]
 value=sum(v[0] for v in vals);err=sum(v[1] for v in vals);lower=float(gammainc(3,x))
 assert value>0 and value>=lower-1e-9
 return dict(j,R=value,quadrature_error=err,positive_lower_R=lower,tail_bound_R=(n+1)*math.exp(-80)*(1+80+80**2/2),scope='floating variance integral; analytic positivity is proved separately')
if __name__=='__main__':
 inp=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True);j=json.loads(inp.read_text())['python'];t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,j))
 (out/'variance72.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-t,rows=rows),indent=2)+'\n');print('PASS',len(rows),'variance integrals')
