from concurrent.futures import ProcessPoolExecutor
from scipy.integrate import quad,IntegrationWarning
import math,json,time,os,warnings

def clip(poly,A,B,C):
 out=[]
 if not poly:return out
 for p,q in zip(poly,poly[1:]+poly[:1]):
  fp=A*p[0]+B*p[1]-C;fq=A*q[0]+B*q[1]-C
  if fp<=0:out.append(p)
  if (fp<0<fq) or (fq<0<fp):
   t=fp/(fp-fq);out.append((p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1])))
 return out

def area(poly):
 return abs(sum(p[0]*q[1]-p[1]*q[0] for p,q in zip(poly,poly[1:]+poly[:1])))/2 if poly else 0

def m0(b,A,r,eps):
 poly=[(0.,0.),(A,0.),(A,r),(0.,r)];t=b*eps/(1-b*eps)
 one=clip(poly,1.,t,t*b);two=clip(poly,t,1.,t*b)
 return area(one)+area(two)-area(clip(one,t,1.,t*b))

def one(a,b,eps,tau):
 lim=min(a,b)
 def integrand(D):
  rr=a-D;h=max(0.,(a*rr-m0(b,a,rr,eps))*eps);h=min(h,b,b-D)
  return 2*((2*b-D)*h-h*h)
 points=[0.,lim]
 for t in [eps,2*eps,eps*eps,2*eps*tau,2*eps*(tau+1),4*eps*(tau+1)]:
  for d in [a-t,b-t]:
   if 0<d<lim:points.append(d)
 points=sorted(set(points));value=error=0.
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter('always',IntegrationWarning)
  for x,y in zip(points,points[1:]):
   v,e=quad(integrand,x,y,epsabs=2e-15,epsrel=1e-12,limit=300)
   value+=v;error+=e
 return value,error,[str(w.message) for w in ws]

def run(pair):
 tau,eps=pair;t=time.monotonic();a=1-eps*tau;b=1+eps*tau
 if a<=0 or b*eps>=1:raise ValueError('invalid geometry')
 f,ef,wf=one(a,b,eps,tau);s,es,ws=one(b,a,eps,tau);val=f+s
 return {'tau':tau,'eps':eps,'u1total':val,'normalized':(val-10*eps/3)/eps**2,'deviation_from_minus_13_over_2':(val-10*eps/3)/eps**2+6.5,'estimated_quad_error_scaled':(ef+es)/eps**2,'warnings':wf+ws,'worker_pid':os.getpid(),'seconds':time.monotonic()-t}
if __name__=='__main__':
 tasks=[(tau,eps) for tau in (0,1,2,5,10,20) for eps in (1e-3,2.5e-4,6.25e-5)]
 with ProcessPoolExecutor(max_workers=6) as ex: rows=list(ex.map(run,tasks))
 print(json.dumps({'method':'scipy quad plus polygon halfplane clip, float64; errors estimated, not certified','scipy_input':tasks,'rows':rows},indent=2))
