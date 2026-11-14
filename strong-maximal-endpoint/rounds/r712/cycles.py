"""Independent deterministic two-cycle expectation integration."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json,math,re
from scipy.integrate import quad
from scipy.special import lambertw
HERE=Path(__file__).resolve().parent
mu=1+float(lambertw(math.exp(-1)).real);r=mu-1
def calc(t):
 def fun(v):
  if v==0:return 1-math.log(abs(t)*(1+r))
  L=(1+r)*(1+t*(1/v-1))
  if abs(L)<1e-6:m=L/2+L*L/6+L**3/12
  elif L==1:m=1
  else:m=1+(1-L)*math.log(abs(1-L))/L
  return -math.log(abs(v))+m
 pts=[-r,0.,1.]
 if t!=1:
  q=-t/(1-t)
  if -r<q<1:pts.append(q)
 q=(1+r)*t/((1+r)*t-r)
 if -r<q<1:pts.append(q)
 pts=sorted(set(pts));pieces=[quad(fun,a,b,epsabs=2e-11,epsrel=2e-11,limit=500,full_output=1) for a,b in zip(pts,pts[1:])]
 return dict(theta=t,mean=sum(p[0] for p in pieces)/(1+r),deficit=2*mu-sum(p[0] for p in pieces)/(1+r),error_estimate=sum(p[1] for p in pieces)/(1+r),quadrature_messages=[str(p[3]) for p in pieces if len(p)>3],subintervals=len(pieces))
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(calc,[s*10.**j for s in (-1,1) for j in (-6,-3,-1,0,1,3,6)]))
 out=dict(workers=12,seed='deterministic',method='conditional exact interval-log integral, scipy QUADPACK, epsabs=epsrel=2e-11',results=rows,scope='finite two-cycle mean diagnostics, not general median bounds')
 (HERE/'cycles_python.json').write_text(json.dumps(out,indent=2)+'\n');print({'cases':len(rows),'all_positive_deficit':all(z['deficit']>0 for z in rows),'messages':sum(len(z['quadrature_messages']) for z in rows)})
