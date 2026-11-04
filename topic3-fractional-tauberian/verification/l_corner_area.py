import numpy as np,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def run(job):
 t,q,N=job;p=2*q-q*q;r=t;ell=np.sqrt(1+(.25+4*p)*r*r)-r;X=(ell+1.5*r)/2;Y=(ell+.5*r)/2;u=1.5*r;v=.5*r;h=q*r
 rects=[(-X+u,X-u,-Y,Y),(-X+h,-X+u,-Y+h,Y-h),(X-u,X-h,-Y+h,Y-h),(-X,-X+h,-Y+v,Y-v),(X-h,X,-Y+v,Y-v)]
 mass=sum((b-a)*(d-c) for a,b,c,d in rects);a0=1-t;nodes,weights=np.polynomial.legendre.leggauss(N);x=(nodes+1)*t;lo=np.zeros_like(x);hi=np.ones_like(x)*2*t;bx=np.unique([z for a,b,c,d in rects for z in (a,b)]);by=np.unique([z for a,b,c,d in rects for z in (c,d)]);sx=2*np.abs(x[None,:]-bx[:,None])
 def val(y):
  L=np.maximum(np.vstack([sx,2*np.abs(y[None,:]-by[:,None])]),1e-30);area=np.zeros_like(L)
  for a,b,c,d in rects:area+=np.maximum(0,np.minimum(b,x+L/2)-np.maximum(a,x-L/2))*np.maximum(0,np.minimum(d,y+L/2)-np.maximum(c,y-L/2))
  return np.max(area/L,axis=0)/np.sqrt(mass)
 for j in range(48):
  md=(lo+hi)/2;ok=val(md)>a0;lo=np.where(ok,md,lo);hi=np.where(ok,hi,md)
 area=4*t*np.dot(weights,(lo+hi)/2)/mass
 return {'t':t,'q':q,'N':N,'area':area,'cubic':(area-2*t*t)/t**3,'mass':mass}
if __name__=='__main__':
 qs=[.1,.2,(4-np.sqrt(7))/6,.3,.4];jobs=[(t,q,1024) for t in [.02,.01,.005,.002,.001] for q in qs];st=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 (Path(__file__).resolve().parents[1]/'results/l_corner_area.json').write_text(json.dumps({'rows':rows,'seconds':time.monotonic()-st,'scope':'unsplit Gaussian boundary quadrature, exploratory only'},indent=2));print([(v['t'],v['q'],v['cubic']) for v in rows]);print(time.monotonic()-st)
