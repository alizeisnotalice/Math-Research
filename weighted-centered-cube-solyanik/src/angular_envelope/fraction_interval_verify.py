#!/usr/bin/env python3
"""15-process Fraction interval certificates for the fixed-d cubic envelope.

For d<2/3 this brackets the unique q_d by signs of h_d, then encloses H(d)
with rational interval arithmetic.  It is a verifier, not a floating solver.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import argparse,hashlib,json,os as _os,platform,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parent;CHECK=ROOT/'checkpoints';VERSION='fraction-envelope-v3-explicit-switch-20261105'

def d_grid():
 ds={Q(i,40) for i in range(1,80)} | {Q(2,3)}
 for k in (3,6,12,24,48,72):
  ds.add(Q(1,2**k));ds.add(Q(2)-Q(1,2**k));ds.add(Q(2,3)-Q(1,2**k));ds.add(Q(2,3)+Q(1,2**k))
 q=Q(3,8);ds.add(3*q*(1-q)*(1-2*q)/(1-3*q+3*q*q))
 return sorted(d for d in ds if 0<d<2)
def r(q):return 1-3*q+3*q*q
def d_of_q(q):return 3*q*(1-q)*(1-2*q)/r(q)
def h(d,q):return r(q)*d-3*q*(1-q)*(1-2*q)
def stationary(d,q):return 6*q**3-3*(3+d)*q*q+3*(1+d)*q-d
def G(d,q):return d*d*(3*q*(1-q)+d*(2*q-1)/2)/(4*q*q*(1-q)*(1-q))
def h_interval(d,left,right):
 # q(1-q) rises on (1/3,1/2); 2q-1 rises there too.
 pL=left*(1-left);pU=right*(1-right)
 nL=3*pL+d*(2*left-1)/2;nU=3*pU+d*(2*right-1)/2
 assert nL>0
 factor=d*d/4
 return factor*nL/(pU*pU),factor*nU/(pL*pL)
def bits_for(d):
 return 128 if min(d,Q(2,3)-d if d<Q(2,3) else Q(1),Q(2)-d)<Q(1,2**12) else 96
def root_certificate(d,bits):
 left,right=Q(1,3),Q(1,2);assert h(d,left)<0<h(d,right)
 for _ in range(bits):
  mid=(left+right)/2
  hm=h(d,mid)
  if hm==0:
   exact=G(d,mid)
   return {'branch':'inner_exact_root','bits':bits,'q_exact':str(mid),'h_exact':'0','stationary_exact':'0','H_exact':str(exact),'valid':exact>0}
  if hm<0:left=mid
  else:right=mid
 lo,hi=h_interval(d,left,right);assert lo<=G(d,(left+right)/2)<=hi
 return {'branch':'inner','bits':bits,'q_lo':str(left),'q_hi':str(right),'q_width':str(right-left),
         'h_lo':str(h(d,left)),'h_hi':str(h(d,right)),'stationary_sign_lo':str(stationary(d,left)),'stationary_sign_hi':str(stationary(d,right)),
         'H_lo':str(lo),'H_hi':str(hi),'H_width':str(hi-lo),'valid':h(d,left)<0<h(d,right) and lo>0}
def outer_certificate(d):
 H=d*(4-d)/(2-d)**2
 return {'branch':'outer','q':str(d/2),'H_exact':str(H),'valid':H>0}
def sig(d):return hashlib.sha256(json.dumps({'d':str(d),'version':VERSION,'rule':'h-sign-bisection+positive-rational-H-interval'},sort_keys=True).encode()).hexdigest()
def task(index,d,barrier=None):
 if barrier is not None:barrier.wait(timeout=30)
 path=CHECK/f'd_{index:03d}.json';s=sig(d)
 if path.exists():
  old=json.loads(path.read_text())
  if old.get('version')==VERSION and old.get('signature')==s:
   old['cache_state']='reused';old['checkpoint_worker_pid']=old.get('worker_pid');old['worker_pid']=_os.getpid();return old
 start=time.time();cert=root_certificate(d,bits_for(d)) if d<Q(2,3) else outer_certificate(d)
 out={'index':index,'d':str(d),'version':VERSION,'signature':s,'worker_pid':_os.getpid(),'checkpoint_worker_pid':_os.getpid(),'cache_state':'computed','seconds':time.time()-start,'certificate':cert,'run_validation':cert['valid']}
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(path);return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);a=ap.parse_args();ds=d_grid();CHECK.mkdir(exist_ok=True)
 from multiprocessing import Manager
 manager=Manager();barrier=manager.Barrier(15);done=[];fail=[];start=time.time();print(json.dumps({'event':'start','d_count':len(ds),'workers':a.workers,'version':VERSION,'platform':platform.platform(),'python':sys.version}),flush=True)
 with ProcessPoolExecutor(max_workers=min(15,a.workers)) as ex:
  # Only the first wave can have 15 concurrent workers; later tasks must not
  # wait on that barrier after workers begin taking a second job.
  fs={ex.submit(task,i,d,barrier if a.workers>=15 and i<15 else None):i for i,d in enumerate(ds)}
  for f in as_completed(fs):
   i=fs[f]
   try:
    z=f.result();done.append(z)
    if not z['run_validation']:fail.append({'index':i,'error':'validation false'})
    print(json.dumps({'event':'done','index':i,'pid':z['worker_pid'],'cache':z['cache_state'],'branch':z['certificate']['branch'],'valid':z['run_validation']}),flush=True)
   except Exception as e:fail.append({'index':i,'error':repr(e)});print(json.dumps({'event':'fail','index':i,'error':repr(e)}),flush=True)
 manager.shutdown();out={'version':VERSION,'d_count':len(ds),'failures':fail,'records':sorted(done,key=lambda x:x['index']),'seconds':time.time()-start,'full_pool':a.workers>=15}
 (ROOT/'fraction_failed.json').write_text(json.dumps(fail,indent=2));(ROOT/'fraction_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps({'event':'finish','d_count':len(done),'failed':len(fail),'seconds':out['seconds']}),flush=True)
 if fail:raise SystemExit(1)
if __name__=='__main__':main()
