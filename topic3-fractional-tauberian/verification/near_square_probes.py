"""Actual finite-contact halo probes; exploratory radial-boundary integration."""
import numpy as np
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json,time,platform
@lru_cache(None)
def quadrature(n):
 z,w=np.polynomial.legendre.leggauss(n);return (z+1)/2,w/2

def value(t,holes,n=128):
 holes=np.asarray(holes,float).reshape(-1,4);mass=1-np.sum((holes[:,1]-holes[:,0])*(holes[:,3]-holes[:,2]));x,w=quadrature(n)
 xe=np.r_[-.5,.5,holes[:,:2].ravel()];ye=np.r_[-.5,.5,holes[:,2:].ravel()]
 def score(r,c):
  X=t*r*x*(-1 if c&1 else 1);Y=t*r*(1-x)*(-1 if c&2 else 1)
  ell=np.concatenate((2*abs(X[None,:]-xe[:,None]),2*abs(Y[None,:]-ye[:,None])),axis=0);ell=np.maximum(ell,1e-15)
  left=X-ell/2;right=X+ell/2;bot=Y-ell/2;top=Y+ell/2
  A=np.maximum(np.minimum(right,.5)-np.maximum(left,-.5),0)*np.maximum(np.minimum(top,.5)-np.maximum(bot,-.5),0)
  for a,b,d,e in holes:
   A-=np.maximum(np.minimum(right,b)-np.maximum(left,a),0)*np.maximum(np.minimum(top,e)-np.maximum(bot,d),0)
  return np.max(A/ell/np.sqrt(mass),axis=0)-(1-t)
 out=0;mins=1.;maxs=0
 for c in range(4):
  lo=np.zeros_like(x);hi=np.full_like(x,2.)
  assert np.min(score(lo,c))>0 and np.max(score(hi,c))<0
  for _ in range(45):
   mid=(lo+hi)/2;inside=score(mid,c)>0;lo=np.where(inside,mid,lo);hi=np.where(inside,hi,mid)
  radius=(lo+hi)/2;out+=np.dot(w,radius**2)/2
  mins=min(mins,float(radius.min()));maxs=max(maxs,float(radius.max()))
 return {'area_ratio':t*t*out/mass,'cubic':(out/mass-2)/t,'radius_min':mins,'radius_max':maxs,'mass_deficit':1-mass}

def holes_model(t,L,amount=1,kind='four'):
 h=amount*t*t/L
 if kind=='edge':return [[-.5,.5,.5-amount*t*t,.5]]
 if kind=='middle':return [[-L/2,L/2,.5-h,.5]]
 rects=[]
 for c in range(4 if kind=='four' else 1):
  a,b=(-.5,-.5+L) if c&1 else (.5-L,.5);d,e=(-.5,-.5+h) if c&2 else (.5-h,.5)
  rects.append([a,b,d,e])
 return rects

def run(args):
 t,L,a,k=args;r=value(t,holes_model(t,L,a,k));return dict(t=t,L=L,amount=a,kind=k,**r)
if __name__=='__main__':
 st=time.monotonic();jobs=[]
 for power in [6,8,10,12]:
  t=2.**-power
  for L in [t,t**.5,.25]:
   for a in [.1,1.,4.]:
    for k in ['four','single','middle']:jobs.append((t,L,a,k))
  jobs.append((t,1.,1.,'edge'))
 with ProcessPoolExecutor(max_workers=6) as p:rows=list(p.map(run,jobs))
 out={'python':platform.python_version(),'numpy':np.__version__,'processes':6,'quadrature_nodes':128,'bisection_steps':45,'randomness':'none; deterministic parameter grid','scope':'finite-contact scale optimization; sampled angular/radial integration, no proof of global or uniform statement','seconds':time.monotonic()-st,'cases':rows};(Path(__file__).resolve().parents[1]/'results/near_square_probes.json').write_text(json.dumps(out,indent=2));print({'cases':len(rows),'seconds':out['seconds'],'largest_cubic':max(r['cubic'] for r in rows)})
