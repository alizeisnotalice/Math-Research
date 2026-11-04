"""Independent eight-scale bounded local exploration; LP is not an infeasibility proof."""
import sys,time,json,platform
from eight_seeded_search import capacity
import numpy as np,scipy
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor

def run(job):
 j,b,t,offset=job;a=1-b;al=1-t;q=(2-np.sqrt(4-12*a*b))/3;gc=(1-4*b*b)/2;gl=2*q-2*q*q;sig=min(gc,gl)+offset;x=1-al*al;r=x/(al+np.sqrt(al*al-4*a*b*x))*(1-sig*t)
 C=r*np.array([[a,b],[-a,b],[a,-b],[-a,-b],[b,a],[-b,a],[b,-a],[-b,-a]])
 rng=np.random.default_rng(20261105+j);bounds=np.array([[-1+1e-9,1/al-1e-9]]*8)
 init=rng.uniform(bounds[:,0],bounds[:,1],size=(160,8))
 models=[np.repeat(2*b,8),np.repeat(-2*b,8),np.r_[np.repeat(-1,4),np.repeat(1-2*q,4)],np.r_[np.repeat(1-2*q,4),np.repeat(-1,4)]]
 # Rationally correct common candidates to next order via exact layer roots.
 for k in range(4):
  model=models[k]
  if k==0:model=np.repeat((2*b*r+np.sqrt(1+4*(a-b)**2*r*r)-1)/t,8)
  if k==1:model=np.repeat((-2*b*r+np.sqrt(1+16*a*b*r*r)-1)/t,8)
  if k>=2:
   p=2*q-q*q;lo=(np.sqrt(1+(1-4*a*b+4*p)*r*r)-r-1)/t;hi=lo+2*(1-q)*r/t
   model=np.r_[np.repeat(lo if k==2 else hi,4),np.repeat(hi if k==2 else lo,4)]
  init[k*32:(k+1)*32]=np.clip(model+rng.normal(size=(32,8))*np.r_[np.repeat(.015,4),np.repeat(.015,4)],bounds[:,0],bounds[:,1]);init[k*32]=np.clip(model,bounds[:,0],bounds[:,1])
 def obj(v):return -capacity(C,1+t*v,al)/t**2
 result=differential_evolution(obj,bounds,init=init,seed=20261105+j,maxiter=420,polish=False,tol=1e-7,atol=1e-8)
 return dict(b=b,t=t,sigma=sig,offset=offset,params=result.x.tolist(),scales=(1+t*result.x).tolist(),slack_over_t2=-result.fun,evaluations=result.nfev,seed=20261105+j,success=bool(result.success),centers=C.tolist())
if __name__=='__main__':
 jobs=[(j,b,.002,off) for j,(b,off) in enumerate(( (b,o) for b in [.1,.2,.25,.265,.35,.45] for o in [-.01,0]))];st=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as p:rows=list(p.map(run,jobs))
 out=dict(rows=rows,seconds=time.monotonic()-st,evaluations=sum(x['evaluations'] for x in rows),workers=8,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,scope='8 independently varying scales, floating search only; no negative certificate')
 open(Path(__file__).resolve().parents[1]/'results/full_direction_search.json','w').write(json.dumps(out,indent=2));print(out)
