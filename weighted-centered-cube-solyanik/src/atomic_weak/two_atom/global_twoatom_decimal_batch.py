#!/usr/bin/env python3
"""15-process high-precision diagnostic for the all-tau two-atom improvement.

Uses u=(r-1)/(r+1) and a logarithmic critical equation, stable as tau->1.
Numerical output does not replace the analytic sign proof.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import hashlib,json,os as _os,platform,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from decimal import Decimal,getcontext
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CHECK=ROOT/'global_checkpoints';VERSION='global-twoatom-decimal-v3-20261105'
TAUS=('1e-12','1e-9','1e-6','1e-4','.001','.01','.05','.1','.25','.5','.75','.9','.99','.999999','0.999999999999')
def P(x,t):return (t*x.ln()).exp()
def K(u,t):
 e=1-t;return -e*u.ln()+(t-(2-t)*u).ln()-((2-t)-t*u).ln()
def solve(t):
 lo=Decimal('1e-100');hi=t/(2-t)*(1-Decimal('1e-100'))
 assert K(lo,t)>0>K(hi,t)
 for _ in range(500):
  m=(lo+hi)/2
  if K(m,t)>0:lo=m
  else:hi=m
 u=(lo+hi)/2;r=(1+u)/(1-u);d=r-1
 J=(P(d,t)+P(2+d,t))/(2*t*(1+d))
 R=(P(1+2*d,t)-P(1-2*d,t)+P(d,t)+P(2-d,t))/(2*t*(1+d))
 # direct F residual at high precision, supplementary to K residual.
 F=r*t*(P(r-1,t-1)+P(r+1,t-1))-(P(r-1,t)+P(r+1,t))
 return u,r,J,R,F,K(u,t)
def source_hash(): return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def sig(i,t):return hashlib.sha256(json.dumps({'i':i,'tau':t,'precision':120,'method':'log-u-bisection500','source_sha256':source_hash(),'python':sys.version,'platform':platform.platform()},sort_keys=True).encode()).hexdigest()
def task(i,ts,barrier):
 barrier.wait(timeout=30);path=CHECK/f'tau_{i:02d}.json';s=sig(i,ts)
 if path.exists():
  old=json.loads(path.read_text())
  if old.get('version')==VERSION and old.get('signature')==s:
   old['cache_state']='reused';old['checkpoint_worker_pid']=old.get('worker_pid');old['worker_pid']=_os.getpid();return old
  stale=CHECK/'stale';stale.mkdir(exist_ok=True)
  path.replace(stale/(path.stem+'_stale_'+str(int(time.time_ns()))+'.json'))
 getcontext().prec=120;t=Decimal(ts);u,r,J,R,F,k=solve(t)
 f43=4*t-1+(4*t-7)*P(Decimal(7),t-1)
 # H is used only on (1/4,1); its displayed derivative is not a sign claim below 1/4.
 if t>Decimal(1)/4:
  hprime=Decimal(7).ln()-Decimal(24)/((7-4*t)*(4*t-1));h_ok=hprime<0;h_value=str(hprime)
 else:
  h_ok=True;h_value=None
 valid=Decimal(1)<r<Decimal(4)/3 and R>J and abs(k)<Decimal('1e-100') and f43<0 and h_ok
 out={'i':i,'tau':ts,'u':str(u),'r':str(r),'tau_J':str(t*J),'tau_R':str(t*R),'tau_gap':str(t*(R-J)),'F_residual':str(F),'log_equation_residual':str(k),'F_4over3_factor':str(f43),'H_prime_on_open_domain':h_value,'worker_pid':_os.getpid(),'checkpoint_worker_pid':_os.getpid(),'cache_state':'computed','version':VERSION,'signature':s,'source_sha256':source_hash(),'python':sys.version,'platform':platform.platform(),'valid':valid}
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(path);return out
def main():
 from multiprocessing import Manager
 CHECK.mkdir(exist_ok=True);manager=Manager();barrier=manager.Barrier(15);done=[];bad=[];print(json.dumps({'event':'start','tasks':15,'version':VERSION,'platform':platform.platform(),'python':sys.version}),flush=True)
 with ProcessPoolExecutor(max_workers=15) as ex:
  fs={ex.submit(task,i,t,barrier):i for i,t in enumerate(TAUS)}
  for f in as_completed(fs):
   i=fs[f]
   try:
    z=f.result();done.append(z)
    if not z['valid']:bad.append({'i':i,'error':'validation'})
    print(json.dumps({'event':'done','i':i,'pid':z['worker_pid'],'cache':z['cache_state'],'valid':z['valid']}),flush=True)
   except Exception as e:bad.append({'i':i,'error':repr(e)});print(json.dumps({'event':'fail','i':i,'error':repr(e)}),flush=True)
 manager.shutdown();out={'version':VERSION,'rows':sorted(done,key=lambda z:z['i']),'failures':bad,'precision':120}
 (ROOT/'global_decimal_failed.json').write_text(json.dumps(bad,indent=2));(ROOT/'global_decimal_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps({'event':'finish','rows':len(done),'failed':len(bad)}),flush=True)
 if bad:raise SystemExit(1)
if __name__=='__main__':main()
