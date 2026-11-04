import numpy as np,json,time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

def run(arg):
 t,k=arg;h=k*t;al=1-t;mass=1-4*h*h
 rects=[(-.5+h,.5-h,-.5,.5),(-.5,-.5+h,-.5+h,.5-h),(.5-h,.5,-.5+h,.5-h)] if h else [(-.5,.5,-.5,.5)]
 nodes,weights=np.polynomial.legendre.leggauss(2048)
 # generous bound; unconditional halo subset [-.5,.5]^2 but radial relevant O(t)
 xmax=3*t; x=(nodes+1)*xmax/2;lo=np.zeros_like(x);hi=np.ones_like(x)*xmax
 def val(y):
  scales=[]
  for a,b,c,d in rects:
   scales.extend([2*np.abs(x-a),2*np.abs(x-b),2*np.abs(y-c),2*np.abs(y-d)])
  L=np.maximum(np.stack(scales),1e-30); area=np.zeros_like(L)
  for a,b,c,d in rects:
   area+=np.maximum(0,np.minimum(b,x+L/2)-np.maximum(a,x-L/2))*np.maximum(0,np.minimum(d,y+L/2)-np.maximum(c,y-L/2))
  return np.max(area/L,axis=0)/np.sqrt(mass)
 for _ in range(52):
  mid=(lo+hi)/2;ok=val(mid)>al;lo=np.where(ok,mid,lo);hi=np.where(ok,hi,mid)
 area=4*xmax/2*np.dot(weights,(lo+hi)/2)/mass
 return [t,k,area,(area-2*t*t)/t**3,area-.5*(1/al-al)**2]
if __name__=='__main__':
 start=time.monotonic();jobs=[(t,k) for t in [.02,.01,.005] for k in [.25,.5,1.,1.5]]
 with ProcessPoolExecutor(max_workers=8) as p:rows=list(p.map(run,jobs))
 out={'rows':rows,'seconds':time.monotonic()-start,'best':[max([r for r in rows if r[0]==t],key=lambda x:x[2]) for t in [.02,.01,.005]],'warning':'unsplit quadrature exploratory only'}
 (Path(__file__).resolve().parents[1]/'results'/'corner_refined.json').write_text(json.dumps(out,indent=2));print(out['best'],out['seconds'])
