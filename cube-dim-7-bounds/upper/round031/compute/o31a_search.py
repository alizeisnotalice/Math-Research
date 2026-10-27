#!/usr/bin/env python3
"""O31A: 12 registered disjoint log-t cells; strict new-model candidate exploration."""
import os,sys,json,time,math,hashlib,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
N=7

def work(i):
 start=time.monotonic();cpu=time.process_time();t_calls=0;rcalls=0;maxrcalls=0;invalid=0;best=None
 low=math.log(2)+(math.log(32)-math.log(2))*i/12;high=math.log(2)+(math.log(32)-math.log(2))*(i+1)/12
 def f(y):
  nonlocal t_calls,rcalls,maxrcalls,invalid,best
  t_calls+=1
  if t_calls>96:raise RuntimeError('96 t evaluations limit')
  t=math.exp(y);h=(1+1/t)/2;A=(1+1/t)**N;B=1-h**N+N*quad(lambda z:1-.5*(1+1/(2*t)-.5*math.exp(-z))**(N-1),0,y,epsabs=1e-12)[0]
  localcalls=0
  def g(R):
   nonlocal localcalls,rcalls
   localcalls+=1;rcalls+=1
   if localcalls>64:raise RuntimeError('64 R calls per t limit')
   s=R/t;D=s**N-((s+1/t)/2)**N
   P=N*quad(lambda z:1-.5*(1+1/(2*t)-.5*math.exp(-z))**(N-1),math.log(R),y,epsabs=1e-12)[0]
   return P-D-math.sqrt(A*(B-D))
  lo=1.;hi=t;gl=g(lo);gh=g(hi)
  if gl<=0 or gh>=0:
   invalid+=1;maxrcalls=max(maxrcalls,localcalls);return 1e100
  for _ in range(50):
   m=(lo+hi)/2
   if g(m)>0:lo=m
   else:hi=m
  R=lo-1e-9
  if R<=1:
   invalid+=1;maxrcalls=max(maxrcalls,localcalls);return 1e100
  gv=g(R);maxrcalls=max(maxrcalls,localcalls)
  if gv<=0:invalid+=1;return 1e100
  sigma=R/t;D=sigma**N-((sigma+1/t)/2)**N;value=(math.sqrt(A)+math.sqrt(B-D))**2
  row={'t_float':repr(t),'R_float':repr(R),'g_float':repr(gv),'F_float':repr(value)}
  if best is None or value<float(best['F_float']):best=row
  return value
 opt=minimize_scalar(f,bounds=(low,high),method='bounded',options={'maxiter':96,'xatol':1e-12})
 if best:
  tq=Q(round(float(best['t_float'])*10**12),10**12);Rq=Q(math.floor(float(best['R_float'])*10**12),10**12)
  best.update({'t_rational':str(tq),'R_rational':str(Rq),'sigma_rational':str(Rq/tq)})
 return {'cell':i,'t_interval':[math.exp(low),math.exp(high)],'pid':os.getpid(),'t_evaluations':t_calls,'R_evaluations_total':rcalls,'R_evaluations_max_per_t':maxrcalls,'no_feasible_evaluations':invalid,'optimizer_success':bool(opt.success),'optimizer_message':str(opt.message),'best':best,'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
if __name__=='__main__':
 start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=12) as ex: rows=list(ex.map(work,range(12)))
 out={'batch':'O31A','method':'12 log-t cells, <=96 evaluations each, <=64 g evaluations per t, 50 bisections, boundary retreat 1e-9','global_optimality_claim':False,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'workers':12,'wall_seconds':time.monotonic()-start,'results':rows};Path(__file__).with_name('o31a_search.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps([{'cell':r['cell'],'best':r['best'],'t_evaluations':r['t_evaluations'],'R_max':r['R_evaluations_max_per_t'],'pid':r['pid']} for r in rows],indent=2))
