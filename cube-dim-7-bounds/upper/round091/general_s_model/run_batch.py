#!/usr/bin/env python3
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import json,math,hashlib,time,resource,multiprocessing as mp,subprocess
from pathlib import Path
from collections import Counter
from scipy.optimize import differential_evolution
P=Path(__file__).resolve().parent
R=json.loads((P/'registration.json').read_text());C=1.1

def build(z):
 s,um,us,uq,v,g=z; mu=s-1+(s+1)*math.log((s+1)/(2*s)); chi=mu+10**(-um); sigma0=chi/(1+chi); sigma=sigma0+10**(-us); q=1-10**(-uq)
 if sigma>=1 or q<=1/s:return None,'domain'
 price=-math.log1p(-sigma); rlo=1.0000001;rhi=s*q-0.0000001
 if rhi<=rlo:return None,'target_band'
 masslog=math.log(rhi/rlo); surplus=masslog-price-C
 if surplus<=0:return None,'nonpositive_cost_surplus'
 gamma=price+C+g*surplus; m=(gamma+masslog)/2; ret=1-min(.001,(m-gamma)/(4*m));e=min(1e-5,(m*ret-gamma)/(30*m),(chi-mu)/20,1/(4*(s+1)))
 if not e>0:return None,'error_budget'
 d=(chi-mu)/4; face=((mu+2*d)/(1-e)+chi)/2;a=chi-face;nu=1/s-(1+1/s)*e;eta=(1-q)/(2*q)
 P0=1-2*e; ulo=P0*(1-e)/(1+e);uhi=(1+e)/(P0*(1-e));M=2*(1+e)*uhi;b=2*ret*(1-e)*ulo-2*e;budget=(1-e)*m*b/M
 low=chi*(1-sigma); delta=low+v*(sigma-low);caprate=nu*eta*(1-sigma)/(2*(1+e));taurate=delta-low
 if min(a,nu,caprate,taurate,sigma-delta,b,budget-gamma,(1-e)*masslog-m)<=0:return None,'secondary_budget'
 return {'s':s,'mu':mu,'chi':chi,'sigma':sigma,'qmax':q,'delta':delta,'gamma':gamma,'price':price,'m':m,'ret':ret,'e':e,'m_face':face,'d_mu':d,'rhi':rhi,'masslog':masslog,'nu':nu,'eta':eta,'M':M,'budget':budget,'caprate':caprate,'taurate':taurate,'slack':gamma-price-C,'a':a},None

def feasible(logL,z):
 L=math.exp(logL);sig=z['sigma'];chi=z['chi'];x=(1-sig)*L;logN=16*logL-math.lgamma(17); logx=math.log(x); e=z['e']; le=math.log(e); d=z['d_mu'];a=z['a'];rhi=z['rhi']; ret=z['ret'];gap=a*x-1
 # All values evaluated in logs where needed, and with conservative N>=2/nu for weighted cap.
 if gap<=0:return False
 mean1=math.log(4*z['s']**4)+2*logx-logN;mean2=math.log(4*z['s']**2)+logx-logN
 mean_max=max(mean1,mean2)+math.log(2)
 checks=[logN>=math.log(400),logN-logx>=math.log(2*rhi),logN-math.log(z['s'])-logx>=math.log(4*(z['s']+1)),
 logN>=math.log(2/z['nu']),2*math.exp(-2*logL)<1e-7,
 math.log(2*(4*z['s']**2+6*z['s']+1))+logx-2*logL<=le,math.log(10*z['s'])-2*logL<=math.log(d),
 math.log(2*z['s']**2)+logx-logN<=le,mean_max<=math.log(e/2),
 math.log(4*16*z['s']**4)+2*logx-logN-2*le<=math.log(e/8),
 math.log(4)-2*le-logN<=math.log(e/8),math.log(z['s']**2)-logN<=math.log(d/2),
 math.log(4*z['s']**4)-2*math.log(d)-logN<=math.log(e/8),
 math.log(2*z['s']**2)+logx-2*math.log(gap)<=math.log(1-ret),
 math.log(z['M'])-z['caprate']*L<=le,
 math.log(rhi)+2*logL-z['taurate']*L<=le,
 math.log(rhi)+2*logL-z['taurate']*L<0,
 -logN<=math.log(z['slack']/4),-z['delta']*L<=math.log(z['slack']/4),-sig*L<=math.log(z['slack']/4),
 logL-(sig-z['delta'])*L<=math.log(z['slack']/4)]
 # target cap (19/20)^(n-1) bounded using n >= max(N,400); use n>=400 only when enough; e may be smaller.
 target_cap_fraction=1-1/(4*rhi)
 min_n=1+math.log(e/8)/math.log(target_cap_fraction)
 checks += [logN>=math.log(min_n),2*math.exp(-2*logL)+math.exp(-min(logN,745))<=3*e/4]
 return all(checks)

def worker(i):
 failures=Counter();calls=0; start=time.monotonic(); cpu0=time.process_time()
 def obj(v):
  nonlocal calls
  calls+=1;z,why=build(v)
  if z is None:failures[why]+=1;return 1e3
  lo=math.log(10**3);hi=math.log(10**4)
  while hi<math.log(2)*200 and not feasible(hi,z):hi+=math.log(2)*4
  if hi>=math.log(2)*200:failures['threshold_above_2^200']+=1;return 999
  for _ in range(42):
   mid=(lo+hi)/2
   if feasible(mid,z):hi=mid
   else:lo=mid
  return hi/math.log(2)
 bounds=[(5+59*i/15,5+59*(i+1)/15),(-.3,6),(1,6),(2,8),(.05,.95),(.05,.95)]
 opt=differential_evolution(obj,bounds,seed=2026102902+i,popsize=7,maxiter=40,polish=False,tol=0,updating='immediate')
 z,why=build(opt.x); usage=resource.getrusage(resource.RUSAGE_SELF)
 out={'worker':i,'pid':os.getpid(),'evaluations':calls,'max_evaluations':1722,'interval':[str(R['intervals'][i]['s_lo']),str(R['intervals'][i]['s_hi'])],'status':'candidate' if z and opt.fun<999 else 'no_candidate','log2_L_model':float(opt.fun),'parameter_decimals':z,'search_coordinates':[float(t) for t in opt.x],'failures':dict(failures),'cpu_seconds':time.process_time()-cpu0,'wall_seconds':time.monotonic()-start,'maxrss_native':usage.ru_maxrss,'maxrss_unit':'bytes on macOS; KiB on Linux','scope':'float diagnostic; must round/register and re-prove all finite gates with rational intervals'}
 (P/f'worker_{i:02}.json').write_text(json.dumps(out,indent=2)+'\n')

def main():
 t=time.monotonic(); workers=[mp.Process(target=worker,args=(i,)) for i in range(15)]
 for w in workers:w.start()
 ids={w.pid for w in workers};peak_cpu=peak_rss=0.;samples=0
 while any(w.is_alive() for w in workers):
  s=subprocess.run(['ps','-Ao','pid=,%cpu=,rss='],capture_output=True,text=True).stdout;cpu=rss=0.
  for row in s.splitlines():
   v=row.split()
   if len(v)==3 and int(v[0]) in ids:cpu+=float(v[1]);rss+=float(v[2])
  peak_cpu=max(peak_cpu,cpu);peak_rss=max(peak_rss,rss);samples+=1;time.sleep(.1)
 for w in workers:w.join()
 rows=[json.loads((P/f'worker_{i:02}.json').read_text()) for i in range(15) if (P/f'worker_{i:02}.json').exists()]; good=sorted([r for r in rows if r['status']=='candidate'],key=lambda r:r['log2_L_model'])
 out={'registration_sha256':hashlib.sha256((P/'registration.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parent_pid':os.getpid(),'worker_pids':sorted(ids),'exit_codes':[w.exitcode for w in workers],'worker_count':len(rows),'total_evaluations':sum(r['evaluations'] for r in rows),'worker_cpu_seconds':sum(r['cpu_seconds'] for r in rows),'wall_seconds':time.monotonic()-t,'sampled_peak_aggregate_cpu_percent':peak_cpu,'sampled_peak_aggregate_rss_KiB':peak_rss,'sample_count':samples,'top_candidates':good[:3],'scope':'general-s float model only, no exact certificate and no Cn theorem'}
 (P/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
