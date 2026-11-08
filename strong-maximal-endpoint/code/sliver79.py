import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq,minimize_scalar
DEN=10**30
def outward(lo,hi):return Q((lo*DEN).__floor__(),DEN),Q((hi*DEN).__ceil__(),DEN)
@lru_cache(None)
def basiclog(t):
 assert 1<=t<=2
 u=(t-1)/(t+1);lo=sum((2*u**(2*k+1)/Q(2*k+1) for k in range(40)),Q(0));hi=lo+2*u**81/(81*(1-u*u))
 return outward(lo,hi)
def logi(t):
 assert t>0
 k=0
 while t<1:t*=2;k-=1
 while t>2:t/=2;k+=1
 lo,hi=basiclog(t);a,b=basiclog(Q(2))
 return outward(lo+k*(a if k>=0 else b),hi+k*(b if k>=0 else a))
def phlog(n,g):
 lo,hi=logi(g);return lo+(n-1)*logi(1+lo)[0],hi+(n-1)*logi(1+hi)[1]
def ff(n,s,r):
 a,b=logi(1+r);c,d=logi(s)
 return (n-1)*r-b-1-d,(n-1)*r-a-1-c
def costbox(n,s,d,lo,hi):
 return logi(2*d/hi)[0]+phlog(n,s*(1+lo))[0],logi(2*d/lo)[1]+phlog(n,s*(1+hi))[1]
def run(z):
 n=z['N'];s=Q(z['sigma']);d=Q(z['d']);rmin=2*d;rmax=2/s-1;assert 0<rmin<=rmax
 fmin=ff(n,s,rmin);fmax=ff(n,s,rmax)
 if fmin[0]>=0:branch='full';rl=rh=rmin
 elif fmax[1]<=0:branch='cap';rl=rh=rmax
 else:
  assert fmin[1]<0 and fmax[0]>0
  branch='interior'
  r0=brentq(lambda r:(n-1)*r-math.log1p(r)-1-math.log(float(s)),float(rmin),float(rmax),xtol=1e-30,rtol=2e-15)
  rq=Q(str(r0));step=max(Q(1,10**24),rq/Q(10**12))
  for _ in range(6):
   rl=max(rmin,rq-step);rh=min(rmax,rq+step)
   if ff(n,s,rl)[1]<0 and ff(n,s,rh)[0]>0:break
   step*=10
  else:raise ArithmeticError('no certified root bracket')
 r=(rl+rh)/2;g=s*(1+r);u=2*d/r;assert 0<u<=1 and s<g<=2 and g*u/(u+2*d)==s
 lo,hi=costbox(n,s,d,rl,rh);grid=phlog(n,s*(1+2*d));cand=costbox(n,s,d,r,r)
 loss=(grid[0]-hi,grid[1]-lo)
 def objective(x):
  v=math.exp(x);lg=math.log(float(s))+math.log1p(v)
  return math.log(float(2*d))-x+lg+(n-1)*math.log1p(lg)
 rstar=float(r);analytic=objective(math.log(rstar))
 direct=minimize_scalar(objective,bounds=(math.log(float(rmin)),math.log(float(rmax))),method='bounded',options={'xatol':1e-12,'maxiter':600})
 assert direct.success and abs(direct.fun-analytic)<2e-6
 asym=None if 'beta' not in z else (z['beta']-1)*math.log(n)-math.log(2*math.e)
 return dict(z,branch=branch,r_lower=str(rl),r_upper=str(rh),r_candidate=str(r),height=str(g),width=str(u),
  log_continuous_optimum_interval=[str(lo),str(hi)],log_candidate_cost_interval=[str(q) for q in cand],log_grid_wholecell_price_interval=[str(q) for q in grid],log_grid_onecolumn_loss_interval=[str(q) for q in loss],
  log_full_grid_fractional_over_continuous_lower=float(logi(s*(1+2*d))[0]-(n-1)*logi(Q(2))[1]-cand[1]),
  grid_sees_no_slab_halo=bool(g*u/(1+d)<1),grid_slab_mean_upper=str(g*u/(1+d)),
  direct_minimizer_success=bool(direct.success),direct_minimizer_objective_gap=direct.fun-analytic,
  asymptotic_log_loss=asym,asymptotic_log_loss_residual=None if asym is None else float((loss[0]+loss[1])/2)-asym)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True)
 ns=[2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576]
 jobs=[dict(N=n,sigma=str(1+Q(1,8*n)),d=str(Q(1,n**b)),beta=b,kind='dimension_scale') for n in ns for b in [2,3,4]]
 jobs +=[dict(N=n,sigma='101/100',d=str(Q(1,10**k)),kind='fixed_dimension',decimal_depth=k) for n in [2,4,8] for k in range(1,9)]
 t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:r=list(pool.map(run,jobs))
 (out/'sliver79.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-t,log_series_terms=40,log_rounding_denominator=str(DEN),rows=r),indent=2)+'\n')
 print('PASS',len(r),'certified continuous optima; max direct gap',max(abs(x['direct_minimizer_objective_gap']) for x in r),'time',time.monotonic()-t)
