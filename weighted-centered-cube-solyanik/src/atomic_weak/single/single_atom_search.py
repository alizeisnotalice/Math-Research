#!/usr/bin/env python3
"""15-process numerical maximization of J_tau(r), r>1 (candidate evidence)."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import hashlib,json,math,os as _os,platform,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CHECK=ROOT/'single_checkpoints';VERSION='singleatom-search-v2-critical-bisection-20261105'
TAUS=(Q(1,512),Q(1,256),Q(1,128),Q(1,64),Q(1,32),Q(1,16),Q(1,12),Q(1,10),Q(1,8),Q(1,7),Q(1,6),Q(1,5),Q(1,4),Q(3,10),Q(1,3))
def J(t,r):return ((r-1)**t+(r+1)**t)/(2*t*r)
def critical(t,r):
 N=(r-1)**t+(r+1)**t;Np=t*((r-1)**(t-1)+(r+1)**(t-1));return r*Np-N
def critical_bisection(t):
 lo=1+1e-14;hi=4.;assert critical(t,lo)>0>critical(t,hi)
 for _ in range(220):
  mid=(lo+hi)/2
  if critical(t,mid)>0:lo=mid
  else:hi=mid
 r=(lo+hi)/2;return r,J(t,r)
def sig(i,t):return hashlib.sha256(json.dumps({'i':i,'tau':str(t),'domain':'r>1','objective':'maximize J','method':'critical_bisection220'},sort_keys=True).encode()).hexdigest()
def task(i,t,barrier):
 barrier.wait(timeout=30);path=CHECK/f'tau_{i:02d}.json';s=sig(i,t)
 if path.exists():
  old=json.loads(path.read_text())
  if old.get('version')==VERSION and old.get('signature')==s:
   old['cache_state']='reused';old['checkpoint_worker_pid']=old.get('worker_pid');old['worker_pid']=_os.getpid();return old
 start=time.time();r,v=critical_bisection(float(t));res=critical(float(t),r)
 out={'i':i,'tau':str(t),'r':r,'J':v,'critical_residual':res,'worker_pid':_os.getpid(),'checkpoint_worker_pid':_os.getpid(),'cache_state':'computed','seconds':time.time()-start,'version':VERSION,'signature':s,'valid':r>1 and math.isfinite(v) and abs(res)<2e-6}
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(path);return out
def main():
 from multiprocessing import Manager
 CHECK.mkdir(exist_ok=True);manager=Manager();barrier=manager.Barrier(15);done=[];bad=[];start=time.time();print(json.dumps({'event':'start','tasks':15,'version':VERSION,'platform':platform.platform(),'python':sys.version}),flush=True)
 with ProcessPoolExecutor(max_workers=15) as ex:
  fs={ex.submit(task,i,t,barrier):i for i,t in enumerate(TAUS)}
  for f in as_completed(fs):
   i=fs[f]
   try:
    z=f.result();done.append(z)
    if not z['valid']:bad.append({'i':i,'error':'invalid'})
    print(json.dumps({'event':'done','i':i,'pid':z['worker_pid'],'cache':z['cache_state'],'valid':z['valid']}),flush=True)
   except Exception as e:bad.append({'i':i,'error':repr(e)});print(json.dumps({'event':'fail','i':i,'error':repr(e)}),flush=True)
 manager.shutdown();out={'version':VERSION,'rows':sorted(done,key=lambda z:z['i']),'failures':bad,'seconds':time.time()-start}
 (ROOT/'single_search_failed.json').write_text(json.dumps(bad,indent=2));(ROOT/'single_search_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps({'event':'finish','rows':len(done),'failed':len(bad),'seconds':out['seconds']}),flush=True)
 if bad:raise SystemExit(1)
if __name__=='__main__':main()
