import numpy as np,json,time,platform
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

def run(job):
 t,d,p,q=job;W=np.exp(d*t);w=1/W;mass=1-4*p*q*t*t;alpha=1-t
 nodes,weights=np.polynomial.legendre.leggauss(192);x=(nodes+1)/2;weights=weights/2
 def val(r):
  xx=t*r*x;yy=t*r*(1-x)
  L=np.maximum(np.stack([2*np.abs(xx-b) for b in [-W/2,W/2,-W/2+p*t,W/2-p*t]]+[2*np.abs(yy-b) for b in [-w/2,w/2,-w/2+q*t,w/2-q*t]]),1e-30)
  bx=np.maximum(0,np.minimum(W/2,xx+L/2)-np.maximum(-W/2,xx-L/2));by=np.maximum(0,np.minimum(w/2,yy+L/2)-np.maximum(-w/2,yy-L/2))
  Ux=sum(np.maximum(0,np.minimum(b,xx+L/2)-np.maximum(a,xx-L/2)) for a,b in [(-W/2,-W/2+p*t),(W/2-p*t,W/2)])
  Uy=sum(np.maximum(0,np.minimum(b,yy+L/2)-np.maximum(a,yy-L/2)) for a,b in [(-w/2,-w/2+q*t),(w/2-q*t,w/2)])
  return np.max((bx*by-Ux*Uy)/L,axis=0)/np.sqrt(mass)
 lo=np.zeros_like(x);hi=np.ones_like(x)*3
 for _ in range(44):
  mid=(lo+hi)/2;ok=val(mid)>alpha;lo=np.where(ok,mid,lo);hi=np.where(ok,hi,mid)
 r=(lo+hi)/2;area=2*t*t*np.dot(weights,r*r)/mass;B=(area-2*t*t)/t**3;R=.5*(1/alpha-alpha)**2
 return [t,d,p,q,area,B,(area-R)/t**3]
if __name__=='__main__':
 start=time.monotonic();jobs=[(t,d,p,q) for t in [.01,.003] for d in [0,.25,.5,.75,.9,.99,1.01] for p in [0,.25,.5,1,1.5] for q in [0,.25,.5,1,1.5]]
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 out={'environment':{'python':platform.python_version(),'numpy':np.__version__},'cases':len(jobs),'rows':rows,'best':[max([r for r in rows if r[0]==t],key=lambda r:r[4]) for t in [.01,.003]],'workers':8,'Gauss_nodes':192,'seconds':time.monotonic()-start,'scope':'finite unsplit quadrature exploration; no upper bound certificate'}
 (Path(__file__).resolve().parents[1]/'results'/'joint_scan.json').write_text(json.dumps(out,indent=2));print(out['best'],out['seconds'])
