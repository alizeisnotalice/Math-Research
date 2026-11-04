#!/usr/bin/env python3
"""120-digit analytic shell-integral diagnostic for the cube single atom.

Numerical confirmation only: it evaluates the supplied exact shell formula;
it neither proves the lower bound nor an upper bound for W_tau.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import hashlib,json,os as _os,platform,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from decimal import Decimal as D,getcontext
from multiprocessing import Manager
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CACHE=ROOT/'cube_checkpoints';VERSION='round27-cube-shell-v2-20261105';PREC=120
TAUS=('0.1','0.01','0.0001','1e-6'); NS=(1,2,3); BS=(1,2,5)
def source_hash():return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def power(x,a):
 if x==0:return D(0)
 assert x>0
 return (a*x.ln()).exp()
def integral_power(lo,hi,e):
 """Integral s^e ds, with positive endpoints except the integrable first shell."""
 q=e+1
 assert q!=0
 return (power(hi,q)-power(lo,q))/q
def shell_mass(n,b,t):
 """Exact integration of s^(t-n) A'(s), breaks t=r-1,r=1+t,r+1."""
 r=1+t;u=r-1 # u=t > 0
 c0=(1+b)*D(2)**(n-1)
 part0=c0*D(n)/t*power(u,t)
 c1=D(2)**(n-1)
 part1=c1*(D(n)*integral_power(u,r,t-1)+b*u*D(n-1)*integral_power(u,r,t-2))
 c2=(D(2)*r)**(n-1)*integral_power(r,r+1,t-n)
 return part0+part1+c2,(part0,part1,c2)
def signature(n,b,ts):
 p={'version':VERSION,'source_sha256':source_hash(),'n':n,'b':b,'tau':ts,'precision':PREC,'python':sys.version,'platform':platform.platform(),'model':'A(s) shell exact derivative'}
 return hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest(),p
def job(item,barrier=None):
 n,b,ts=item;path=CACHE/f'n{n}_b{b}_tau{ts.replace(".","p")}.json';sig,payload=signature(n,b,ts)
 if barrier is not None:barrier.wait(timeout=60)
 if path.exists():
  old=json.loads(path.read_text())
  if old.get('signature')==sig:
   old['cache_state']='reused';old['replay_pid']=_os.getpid();return old
  stale=CACHE/'stale';stale.mkdir(exist_ok=True);path.replace(stale/(path.stem+'_stale_'+str(time.time_ns())+'.json'))
 getcontext().prec=PREC;t=D(ts);bb=D(b);mu,parts=shell_mass(n,bb,t);r=1+t;V=D(2)**n;cost=V*r**n;R=mu/cost;target=D(n)*(1+bb)/2
 # Independent comparisons requested by the three-scale framework.
 ea=(1+bb)/2;ball_lower=D(n)*ea/t*power(t,t)/(1+t)**n
 eta=power(t,D(1)/(n+1));bk=D(2)**n
 tri_upper=bk*bb*power(D(3),D(n)-t)+D(n)*(D(2)**n)*ea/(D(n)-t)*power(1+1/eta,D(n)-t)+D(n)*ea/t*power(eta,t)*power(1+eta,D(n)-t)
 valid=mu>0 and R>0 and R>=ball_lower and R<=tri_upper
 out={'n':n,'b':b,'tau':ts,'r':str(r),'shell_mass':str(mu),'cost':str(cost),'R':str(R),'tauR':str(t*R),'target':str(target),'target_error':str(t*R-target),'ball_lower':str(ball_lower),'three_scale_eta':str(eta),'three_scale_upper':str(tri_upper),'parts':[str(x) for x in parts], 'worker_pid':_os.getpid(),'cache_state':'computed','signature':sig,'signature_payload':payload,'valid':valid}
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(path);return out
def main():
 CACHE.mkdir(exist_ok=True);items=[(n,b,t) for n in NS for b in BS for t in TAUS];rows=[];bad=[]
 print(json.dumps({'event':'start','tasks':len(items),'workers':15,'precision':PREC,'taus':TAUS}),flush=True)
 with Manager() as manager:
  # First wave has one task per worker and a barrier: this records 15 actual PIDs.
  barrier=manager.Barrier(15)
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={ex.submit(job,x,barrier):x for x in items[:15]}
   for f in as_completed(fs):
    try:
     z=f.result();rows.append(z);print(json.dumps({'event':'done','n':z['n'],'b':z['b'],'tau':z['tau'],'pid':z['worker_pid'],'tauR':z['tauR'],'cache':z['cache_state']}),flush=True)
    except Exception as e:bad.append({'task':fs[f],'error':repr(e)});print(json.dumps({'event':'failure','task':fs[f],'error':repr(e)}),flush=True)
   fs={ex.submit(job,x):x for x in items[15:]}
   for f in as_completed(fs):
    try:
     z=f.result();rows.append(z);print(json.dumps({'event':'done','n':z['n'],'b':z['b'],'tau':z['tau'],'pid':z['worker_pid'],'tauR':z['tauR'],'cache':z['cache_state']}),flush=True)
    except Exception as e:bad.append({'task':fs[f],'error':repr(e)});print(json.dumps({'event':'failure','task':fs[f],'error':repr(e)}),flush=True)
 rows.sort(key=lambda z:(z['n'],z['b'],D(z['tau'])));out={'version':VERSION,'precision':PREC,'rows':rows,'failures':bad,'unique_inputs':len(items),'configured_workers':15,'actual_worker_pids':sorted(set(z['worker_pid'] for z in rows))}
 (ROOT/'cube_shell_summary.json').write_text(json.dumps(out,indent=2));(ROOT/'cube_shell_failures.json').write_text(json.dumps(bad,indent=2));print(json.dumps({'event':'finish','rows':len(rows),'failures':len(bad)}),flush=True)
 if bad:raise SystemExit(1)
if __name__=='__main__':main()
