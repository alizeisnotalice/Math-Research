#!/usr/bin/env python3
"""C34A: adaptive candidate search and exact rational lower certificates."""
import os,sys,json,time,math,hashlib,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
from scipy.optimize import differential_evolution
YS=(14,16,20,24,28,32,40,48,64,96,128,256)
def exact(Y,a,b,c):
 Y=Q(Y); a=Q(a);b=Q(b);c=Q(c)
 assert Q(1,10)*b<=a<=Q(9,10)*b and Q(2,25)<=b<=Q(4,5) and Q(1,200)<=c<=Q(2,5)
 assert (b+c)*Y>=1
 f1=1-4*(b+c)/Y; f2=1-2*(b+c)**2/(1-b/Y);assert f1>0 and f2>0
 KA=c/(4*(1+c/Y));CL=f1*f2;r=Q(9,32)*a*Y
 term=Q(1);S=term
 for k in range(1,129):term=term*r/k;S+=term
 # e^r > S for r>0, so 1/S is a strict upper bound on exp(-r).
 KP=(b-a)*(CL-Q(1,2)/S)-1/Y-Q(3,8*Y**3)
 return KA,KP,f1,f2

def baseline():
 a,b,c=Q(1,10),Q(29,100),Q(21,100);KA,KP,f1,f2=exact(64,a,b,c);q=Q(1,20)
 return {'Y':64,'a':str(a),'b':str(b),'c':str(c),'q':str(q),'KA_approx':float(KA),'KP_certified_lower_approx':float(KP),'q_lt_both':q<KA and q<KP,'rational_strict_KA_margin':str(KA-q),'rational_strict_KP_margin':str(KP-q)}
def run(Y):
 start=time.monotonic();cpu=time.process_time();calls=0; valid=0
 def fun(v):
  nonlocal calls,valid
  calls+=1
  if calls>20000:raise RuntimeError('objective budget exceeded')
  theta,b,c=v;a=theta*b
  if (b+c)*Y<1 or b>=Y:return 1e3
  f1=1-4*(b+c)/Y;f2=1-2*(b+c)**2/(1-b/Y)
  if f1<=0 or f2<=0:return 1e3
  KA=c/(4*(1+c/Y));KP=(b-a)*(f1*f2-.5*math.exp(-9*a*Y/32))-1/Y-3/(8*Y**3)
  if KP<=0:return 1e3
  valid+=1;return -min(KA,KP)
 result=differential_evolution(fun,[(.1,.9),(.08,.8),(.005,.4)],seed=349001+Y,popsize=12,maxiter=500,tol=1e-9,atol=0,polish=False,workers=1,updating='immediate')
 theta,b0,c0=map(float,result.x);a0=theta*b0;status='certified'
 try:
  a,b,c=(Q(round(v*10**6),10**6) for v in (a0,b0,c0));KA,KP,f1,f2=exact(Y,a,b,c)
  assert KP>0; bound=min(KA,KP);q=Q((bound.numerator*10**12//bound.denominator)-1,10**12)
  assert 0<q<KA and q<KP
  cert={'a':str(a),'b':str(b),'c':str(c),'q':str(q),'KA_lower':str(KA),'KP_lower':str(KP),'strict_KA_margin':str(KA-q),'strict_KP_margin':str(KP-q),'KA_approx':float(KA),'KP_lower_approx':float(KP),'q_approx':float(q),'factor1_positive':True,'factor2_positive':True,'exponential_upper':'1/sum_{k=0}^{128} r^k/k!, r=9aY/32'}
 except (AssertionError,ValueError,ZeroDivisionError) as e:status='no certified candidate';cert={'error':repr(e)}
 return {'Y':Y,'pid':os.getpid(),'seed':349001+Y,'raw_theta':repr(theta),'raw_a':repr(a0),'raw_b':repr(b0),'raw_c':repr(c0),'raw_objective':repr(-float(result.fun)),'status':status,'certificate':cert,'objective_calls':calls,'valid_calls':valid,'max_calls':20000,'stop_message':str(result.message),'iterations':int(result.nit),'success_flag':bool(result.success),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'candidate only; no global optimality claim; conditional on independent theorem gate L>=max(64,Y^2)'}
if __name__=='__main__':
 path=Path(__file__);base=baseline();print('BASELINE',json.dumps({k:v for k,v in base.items() if 'margin' not in k}),flush=True)
 if '--baseline-only' in sys.argv:exit()
 start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=12) as ex:results=list(ex.map(run,YS))
 out={'batch':'C34A','model_status':'candidate constants only, not an operator bound before theorem ACCEPT','L_requirement':'L >= max(64,Y^2)','script_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'workers_max':12,'wall_seconds':time.monotonic()-start,'baseline':base,'results':results}
 path.with_name('c34a_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([{k:r[k] for k in ('Y','status','objective_calls','raw_objective','pid','wall_seconds')} for r in results],indent=2))
