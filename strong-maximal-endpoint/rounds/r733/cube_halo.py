#!/usr/bin/env python3
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[name]='1'
import json,time,math,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import mpmath as mp
import scipy
from scipy import special,integrate
HERE=Path(__file__).resolve().parent
DIMS=[1,2,4,8,16,32,64,128,256,1024,4096,65536]
def calc(task):
 n,dps=task;mp.mp.dps=dps+20;ll=mp.log(2);mode=mp.sqrt(2*n*ll);cut=min(n,int(mp.ceil(3*mode+100)));coef=mp.mpf(1);terms=[]
 for j in range(1,cut+1):
  coef*=((n-j+1)*2*ll)/(j*j);terms.append(coef*mp.hyp1f1(j,j+1,ll))
 g=1+mp.fsum(terms)
 if cut<n:
  nextcoef=coef*(n-cut)*2*ll/(cut+1)**2;ratio=2*ll*(n-cut-1)/((cut+1)*(cut+2));tail=2*nextcoef/(1-ratio)
 else:tail=mp.mpf(0)
 out=dict(n=n,target_digits=dps,guard_digits=20,positive_terms_used=cut,log_halo=mp.nstr(mp.log(g),dps),relative_analytic_tail_bound_float=mp.nstr(tail/g,dps),log_halo_divided_by_n=mp.nstr(mp.log(g)/n,20),coarse_log_upper=mp.nstr(ll+2*mode,20))
 if dps==40:
  val,err=integrate.quad(lambda s:2*math.exp(s)*special.eval_genlaguerre(n-1,1,-2*s),0,math.log(2),epsabs=1e-9,epsrel=1e-11,limit=200)
  out.update(independent_scipy_log_halo=math.log1p(val),scipy_integral_relative_error_estimate=err/(1+val),scipy_minus_mp_log=math.log1p(val)-float(mp.log(g)))
 return out
if __name__=='__main__':
 start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:r=list(pool.map(calc,[(n,p) for n in DIMS for p in [40,70]]))
 out=dict(status='complete',scope='Cube full strong-maximal half-level halo only, not generalCn; finite floating checks of exact positive formula',workers=12,precisions=[40,70],guard_digits=20,dimensions=DIMS,scipy_version=scipy.__version__,mpmath_version=mp.__version__,python=sys.version.split()[0],wall_seconds=time.monotonic()-start,results=r)
 (HERE/'cube_halo.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status='complete',cases=len(r),wall_seconds=out['wall_seconds'],maximum_scipy_log_difference=max(abs(x.get('scipy_minus_mp_log',0)) for x in r),examples=[[x['n'],x['log_halo']] for x in r if x['target_digits']==70 and x['n'] in [1,16,256,4096,65536]])))
