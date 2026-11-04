"""Four independent rectangular corner holes; numerical cubic functional."""
import numpy as np,json,time,platform,scipy
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from functools import lru_cache

@lru_cache(maxsize=4)
def quadrature(N):
 nodes,weights=np.polynomial.legendre.leggauss(N)
 return (nodes+1)/2,weights/2


def value(v,N=128):
 k=np.asarray(v).reshape(4,2);mass=k[:,0]*k[:,1];M=mass.sum();x,weights=quadrature(N);m=np.minimum(x,1-x);result=0
 for c in range(4):
  ix,iy,ij=c^1,c^2,c^3
  candidates=np.vstack([-m,m,x-k[ix,0],1-x-k[iy,1],x-k[ij,0],1-x-k[ij,1]])
  valid=(candidates>=-m-1e-14)&(candidates<=m+1e-14);u=np.clip(candidates,-m,m);s=x-u;q=1-x-u
  G=u*u+u+x*(1-x)+M/2-mass[c]-np.maximum(k[ix,0]-s,0)*k[ix,1]-k[iy,0]*np.maximum(k[iy,1]-q,0)-np.maximum(k[ij,0]-s,0)*np.maximum(k[ij,1]-q,0)
  result+=np.dot(weights,np.max(np.where(valid,G,-np.inf),axis=0))
 return float(result)

def run(j):
 rng=np.random.default_rng(20261105+j);init=rng.uniform(0,1.5,size=(80,8));init[0]=0;init[1]=1;init[2]=np.array([1,1,0,0,0,0,0,0]);init[3]=np.array([1,1,1,1,0,0,0,0])
 opt=differential_evolution(lambda v:-value(v),[(0,2)]*8,init=init,seed=20261105+j,maxiter=220,polish=False,tol=1e-7,atol=1e-8)
 return dict(seed=20261105+j,params=opt.x.tolist(),B=-opt.fun,evaluations=opt.nfev,success=bool(opt.success),B1024=value(opt.x,1024))
if __name__=='__main__':
 # Exploratory integration is unsplit at x=1/2; exact integrator supplies regressions.
 st=time.monotonic()
 with ProcessPoolExecutor(max_workers=6) as p:rows=list(p.map(run,range(6)))
 out=dict(rows=rows,seconds=time.monotonic()-st,evaluations=sum(r['evaluations'] for r in rows),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,scope='four independent rectangular holes; unsplit Gaussian cubic integration, exploratory only');(Path(__file__).resolve().parents[1]/'results/general_corners_search.json').write_text(json.dumps(out,indent=2));print({'evaluations':out['evaluations'],'seconds':out['seconds'],'best_B1024':max(r['B1024'] for r in rows)})
