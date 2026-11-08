import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,itertools,json,time,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq,differential_evolution
C=1.2784645427610738

def halo(params,mode='finite',mixed=False):
 xy=np.array([[0.,0.],params[:2],params[2:4]],float);u=xy.sum(1);z=np.array([0.,params[4],params[5]]);w=np.exp(z-z.max());w/=w.sum();center=float(w@u)/2;xy-=center;u=xy.sum(1);Q=math.exp(params[6]);c=C if mode=='finite' else 0.;labels=[0,1,2] if mixed else [0,0,0]
 def H(s):return s+2*max(s-c,0)+max(s-2*c,0) if mode=='finite' else s
 logt=-Q-math.log(H(Q));roots=[]
 for bits in range(1,8):
  ids=[i for i in range(3) if bits>>i&1];lo=-float(min(u[ids]))+1e-9;hi=Q-float(min(u[ids]))+20
  def f(t):return math.log(sum(w[i]/H(t+u[i]) for i in ids))-t-logt
  roots.append(brentq(f,lo,hi,xtol=3e-10,rtol=1e-12))
 low=-float(max(u));high=max(roots);cuts=[low,high]+roots
 for sig in (range(4) if mode=='finite' else [0]):
  js=[sig^t for t in labels]
  for i in range(3):
   for k in range(3):cuts.append(c*((js[i]>>1&1)+(js[k]&1))-xy[i,1]-xy[k,0])
 cuts=sorted(set(t for t in cuts if low<=t<=high));total=0.
 for l,r in zip(cuts,cuts[1:]):
  if r-l<1e-11:continue
  mid=(l+r)/2;al=0;be=0.
  for sig in (range(4) if mode=='finite' else [0]):
   events=[]
   for i in range(3):
    J=sig^labels[i];left=-mid-xy[i,1]+c*(J>>1&1);right=xy[i,0]-c*(J&1)
    if left<right:
     amp=w[i]*math.exp(-mid)/H(mid+u[i]);events.extend([(left,amp,-1.,-xy[i,1]+c*(J>>1&1)),(right,-amp,0.,right)])
   events.sort();cur=0.
   for e,f in zip(events,events[1:]):
    cur+=e[1]
    if cur>math.exp(logt)*(1+2e-13):al+=f[2]-e[2];be+=f[3]-e[3]
  h=r-l;em=math.expm1(h);second=h+(h-1)*em if h>=1e-4 else h*h/2+h**3/3+h**4/8+h**5/30
  total+=math.exp(l-Q)/H(Q)*((al*l+be)*em+al*second)
 return float(total)

def run(job):
 R=job['R'];start=time.monotonic();mode=job['mode'];mixed=job['mixed'];seed=job['seed'];best=[-1,None];invalid=[0]
 def obj(x):
  try:val=halo(x,mode,mixed)
  except (ValueError,OverflowError,ZeroDivisionError):invalid[0]+=1;return 1e3
  if val>best[0]:best[:]=[val,x.tolist()]
  return -val
 def stop(x,conv):return time.monotonic()-start>25
 bounds=[(-R,R)]*4+[(-4,4)]*2+[(math.log(.05),math.log(32))]
 res=differential_evolution(obj,bounds,seed=seed,popsize=6,maxiter=160,tol=1e-6,polish=False,callback=stop,workers=1)
 return {'spec':job,'ratio':best[0],'parameters':best[1],'evaluations':res.nfev,'invalid_objective_evaluations':invalid[0],'solver_success':bool(res.success),'message':str(res.message),'seconds':time.monotonic()-start,'scope':'floating three-atomic-source diagnostic, no directed root or halo certificate'}
if __name__=='__main__':
 jobs=[{'R':R,'mode':md,'mixed':mix,'seed':7100+i*2+k} for i,(R,md,mix) in enumerate([(4,'finite',False),(4,'finite',True),(12,'finite',False),(12,'finite',True),(4,'one_sided',False),(12,'one_sided',False)]) for k in [0,1]]
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'points71.json').write_text(json.dumps({'seconds':time.monotonic()-t,'workers':12,'rows':rows},indent=2)+'\n');print('Completed 12 point-source optimizations; maximum',max(r['ratio'] for r in rows))
