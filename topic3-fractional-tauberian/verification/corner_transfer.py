import numpy as np,json,time
from scipy.integrate import quad
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

def G(x,k):
 y=1-x;m=min(x,y)
 aa=[-2*m,2*m,np.clip(2*(x-k),-2*m,2*m),np.clip(2*(y-k),-2*m,2*m)]
 return max(a*a/4+a/2+x*y+2*k*k-(k+max(0,k-x+a/2))*(k+max(0,k-y+a/2)) for a in aa)
def B(k):
 if k<=.5:return 2-4*k**3
 if k<=2/3:return 2+4*(k-1)*k*k
 if k<=1:return 6-18*k+22*k*k-8*k**3
 return -2+12*k-8*k*k
def run(arg):
 t,k=arg;t=np.longdouble(t);k=np.longdouble(k);h=k*t;mass=1-4*h*h
 R=np.array([[-.5+h,.5-h,-.5,.5],[-.5,-.5+h,-.5+h,.5-h],[.5-h,.5,-.5+h,.5-h]],dtype=np.longdouble)
 def q(x):
  x=np.longdouble(x);lo=np.longdouble(-50);hi=np.longdouble(50)
  for _ in range(62):
   mid=(lo+hi)/2;r=1+t*mid;z=t*r*np.array([x,1-x]);L=np.unique(np.concatenate([2*abs(z[0]-R[:,0]),2*abs(z[0]-R[:,1]),2*abs(z[1]-R[:,2]),2*abs(z[1]-R[:,3])]))
   ar=np.maximum(0,np.minimum(R[:,1,None],z[0]+L/2)-np.maximum(R[:,0,None],z[0]-L/2))*np.maximum(0,np.minimum(R[:,3,None],z[1]+L/2)-np.maximum(R[:,2,None],z[1]-L/2))
   v=np.max(ar.sum(axis=0)/L)/np.sqrt(mass)
   if v>1-t:lo=mid
   else:hi=mid
  return float((lo+hi)/2)
 # integrate q-G, not a tiny area difference; residual is expected to vanish
 residue,error=quad(lambda x:q(x)-G(x,float(k)),0,.5,epsabs=2e-7,limit=160,points=[.1,.2,.3,.4])
 return {'t':float(t),'k':float(k),'B':B(float(k)),'scaled_radial_residual':8*residue,'quadrature_error_estimate':8*error}
if __name__=='__main__':
 start=time.monotonic();jobs=[(t,k) for k in [.25,.5,.625,.75,1.,1.5] for t in [.01,.003,.001]]
 with ProcessPoolExecutor(max_workers=8) as p:rows=list(p.map(run,jobs))
 out={'rows':rows,'seconds':time.monotonic()-start,'workers':8,'dtype':'numpy.longdouble','epsilon':float(np.finfo(np.longdouble).eps),'scope':'finite floating point convergence evidence, adaptive error estimates not certificates'}
 (Path(__file__).resolve().parents[1]/'results'/'corner_transfer_python.json').write_text(json.dumps(out,indent=2));print(out)
