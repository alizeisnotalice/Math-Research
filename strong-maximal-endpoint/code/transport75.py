import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.linalg import cholesky
import numpy as np
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
def params(n):
 rho=math.sqrt(math.log(n)/n);cr=MU-math.sqrt(MU)/2*rho;lr=math.log1p(math.exp(-cr));return cr,cr-rho,cr+rho,math.log1p(math.exp(-(cr-rho)))-lr,math.log1p(math.exp(-(cr+rho)))-lr

def run(j):
 n=j['N'];kind=j['kind'];ref=j['refinement'];q={'half':.5,'sixteenth':1/16,'inverseN':1/n,'inverseN2':1/n**2}.get(kind,0.);cr,cl,ch,lm,lp=params(n);lo=-lm/(1-q);hi=-lp/(1-q)
 if q>0:
  eps=q/n
  def f(u):return math.exp(n*q*u-.5*(n*u)**2)
  points=[v for v in [0,q/n,-8/n,8/n] if lo<v<hi];rec,er=quad(f,lo,hi,points=sorted(set(points)),epsabs=1e-13*ref/n,epsrel=2e-11*ref,limit=160)
  vlo=q*lo;vhi=q*hi
  def g(v):return math.exp(n*v-.5*(v/eps)**2)
  points=[v for v in [0,n*eps**2,-8*eps,8*eps] if vlo<v<vhi];img,ei=quad(g,vlo,vhi,points=sorted(set(points)),epsabs=1e-13*ref*q/n,epsrel=2e-11*ref,limit=160)
  ratio=rec/img;assert abs(q*ratio-1)<1e-9
  return dict(j,q=q,c_reference=cr,c_lower=cl,c_upper=ch,u_lower=lo,u_upper=hi,receiver_weighted_integral=rec,image_integral=img,pushforward_ratio=ratio,predicted_ratio=1/q,receiver_quad_error_estimate=er,image_quad_error_estimate=ei)
 eps=1/(n*2**j['tube_power']);receiver=hi-lo;target,err=quad(lambda v:math.exp(n*v-.5*(v/eps)**2),-12*eps,12*eps,epsabs=1e-13*ref*eps,epsrel=2e-11*ref,points=[0],limit=160);exact=eps*math.sqrt(2*math.pi)*math.exp((n*eps)**2/2)
 assert abs(target/exact-1)<1e-10
 return dict(j,q=0.,c_reference=cr,c_lower=cl,c_upper=ch,u_lower=lo,u_upper=hi,receiver_weighted_integral=receiver,image_integral=target,pushforward_ratio=receiver/target,Gaussian_full_image_integral=exact,image_quad_error_estimate=err)

def matrix(j):
 n=j['N'];kind=j['kind'];q={'half':.5,'sixteenth':1/16,'inverseN':1/n,'inverseN2':1/n**2}.get(kind,0.);x=np.arange(1,n+1,dtype=float);B=np.eye(n)-(1-q)*np.outer(x,x)/(x@x);ev=np.linalg.eigvalsh(B);ld=2*float(np.log(np.diag(cholesky(B,lower=True))).sum()) if q>0 else None;sgn=1. if q>0 else 0.;err=abs(ev[0]-q);assert err<1e-11 and abs(ev[-1]-1)<1e-11
 return dict(j,q=q,minimum_eigenvalue=float(ev[0]),maximum_eigenvalue=float(ev[-1]),det_sign=float(sgn),log_abs_det=ld,det_log_error=None if not q else float(ld-math.log(q)))
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);ns=[2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576];jobs=[dict(N=n,kind=k,refinement=r) for n in ns for k in ['half','sixteenth','inverseN','inverseN2'] for r in [1,.1]]+[dict(N=n,kind='collapse',tube_power=t,refinement=r) for n in ns for t in [0,4,8,16] for r in [1,.1]];mj=[dict(N=n,kind=k) for n in ns[:8] for k in ['half','sixteenth','inverseN','inverseN2','collapse']];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs));mr=list(pool.map(matrix,mj))
 (out/'transport75.json').write_text(json.dumps(dict(mu=MU,workers=12,seconds=time.monotonic()-start,rows=rows,matrices=mr,scope='radial reanchoring geometry only; probes are not actual qualified source-cost fields'),indent=2)+'\n');print('PASS',len(rows),'transport integrals and',len(mr),'Jacobian matrix checks')
