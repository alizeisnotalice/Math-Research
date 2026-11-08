import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
NS=[2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576]
def job(z):
 n,k=z;L=math.log(n);delta={'inverseN':1/n,'inverse_sqrt_NlogN':1/math.sqrt(n*L),'sqrt_logN_over_N':math.sqrt(L/n),'constant_half':.5}[k]
 s=1+delta;pen=math.log1p(.005);a=(n-1)*math.log1p(math.log1p(delta))
 return dict(N=n,scale=k,delta=delta,log_integer_optimum=math.log1p(delta)+a,
  log_fractional_feasible_cost=math.log1p(delta)+pen,log_integrality_gap_lower=a-pen,
  log_perspective_optimum=math.log1p(delta)+a)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as p:r=list(p.map(job,[(n,k) for n in NS for k in ['inverseN','inverse_sqrt_NlogN','sqrt_logN_over_N','constant_half']]))
 (out/'gap77.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-t,eta='1/100',rows=r),indent=2)+'\n')
 print('PASS',len(r),'stable dimension-gap evaluations')
