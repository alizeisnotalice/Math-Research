#!/usr/bin/env python3
"""Float-only 15-process counterexample search for finite rational vectors.

The verifier is deliberately separate. Inputs generated here are positive
rational (represented by Fraction before conversion), have exact mean one,
and include coarse-grained near-critical families.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[_k]='1'
import argparse, hashlib, json, os as _os, platform, sys, time
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parent;CHECK=ROOT/'checkpoints';VERSION='multi-p3-search-v4-law-dedupe-20261105'
DIMS=(3,4,8,16,64); C_FLOAT=(11+5*5**0.5)/8

def normalize(w,x):
 m=sum((a*b for a,b in zip(w,x)),Q(0));return w,[z/m for z in x]
def random_vector(n,seed):
 raww=[Q((11*seed+7*i)%97+1,1) for i in range(n)];w=[z/sum(raww,Q(0)) for z in raww]
 rawx=[Q((17*seed+13*i)%89+1,(seed+3*i)%23+2) for i in range(n)]
 return normalize(w,rawx)
def critical_vector(n,q,a):
 lo=max(1,n//2);hi=n-lo;w=[q/Q(lo) for _ in range(lo)]+[(1-q)/Q(hi) for _ in range(hi)]
 b=(1-q*a)/(1-q); x=[]
 for i in range(lo): x.append(a*(Q(9,8) if i%2==0 else Q(7,8)))
 # Correct the final low point to make its group mean exactly a.
 x[-1]=lo*a-sum(x[:-1],Q(0))
 y=[]
 for i in range(hi): y.append(b*(Q(10,9) if i%2==0 else Q(8,9)))
 y[-1]=hi*b-sum(y[:-1],Q(0));return w,x+y
def embedded_binary_vector(n,q,a):
 lo=max(1,n//2);hi=n-lo;w=[q/Q(lo) for _ in range(lo)]+[(1-q)/Q(hi) for _ in range(hi)]
 b=(1-q*a)/(1-q);return w,[a for _ in range(lo)]+[b for _ in range(hi)]
def equal_atom_vector(n):
 w=[Q(1,n) for _ in range(n)];x=[Q(1) for _ in range(n)];x[0]=Q(1,2);x[1]=Q(3,2);return w,x
def law_signature(w,x):
 merged={}
 for r,z in zip(w,x):merged[z]=merged.get(z,Q(0))+r
 return '|'.join('%s:%s' % (z,merged[z]) for z in sorted(merged))
def eval_float(w,x):
 t=sum(float(a)*float(z)**3 for a,z in zip(w,x))-1.0;d=sum(float(a)*abs(float(z)-1.0) for a,z in zip(w,x));return t/(d*d),t,d
def sig(job):return hashlib.sha256(json.dumps({'job':job,'dims':DIMS,'randomSeeds':12,'critical':'embedded+spread,q381966±2e-6,a2^-8,-16,-32;equalatom','law':'merge equal f then sort(value,totalmass)'},sort_keys=True).encode()).hexdigest()
def task(job,barrier=None):
 if barrier is not None:barrier.wait(timeout=30)
 path=CHECK/f'job_{job:02d}.json';s=sig(job)
 if path.exists():
  old=json.loads(path.read_text())
  if old.get('version')==VERSION and old.get('signature')==s:
   old['cache_state']='reused';old['checkpoint_worker_pid']=old.get('worker_pid');old['worker_pid']=_os.getpid();return old
 best=None;bad=[];count=0;laws=set();start=time.time()
 for n in DIMS:
  candidates=[]
  candidates += [('random',*random_vector(n,1000+n*31+k)) for k in range(12)]
  for q in (Q(381964,10**6),Q(381966,10**6),Q(381968,10**6)):
   for k in (8,16,32):
    candidates.append(('critical_embedded',*embedded_binary_vector(n,q,Q(1,2**k))))
    candidates.append(('critical_spread',*critical_vector(n,q,Q(1,2**k))))
  candidates.append(('equal_atom',*equal_atom_vector(n)))
  for index,(kind,w,x) in enumerate(candidates):
   if index%15!=job:continue
   laws.add(law_signature(w,x))
   value,t,d=eval_float(w,x);count+=1
   if not (d>0 and t>0 and value>=C_FLOAT-1e-10):bad.append({'n':n,'kind':kind,'value':value})
   if best is None or value<best['value']:best={'n':n,'kind':kind,'value':value,'T':t,'d':d}
 out={'job':job,'version':VERSION,'signature':s,'worker_pid':_os.getpid(),'checkpoint_worker_pid':_os.getpid(),'cache_state':'computed','seconds':time.time()-start,'evaluations':count,'law_signatures':sorted(laws),'best':best,'bad':bad,'run_validation':not bad}
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(path);return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);a=ap.parse_args();CHECK.mkdir(exist_ok=True)
 from multiprocessing import Manager
 manager=Manager();barrier=manager.Barrier(15);done=[];failed=[];start=time.time();print(json.dumps({'event':'start','tasks':15,'workers':a.workers,'version':VERSION,'platform':platform.platform(),'python':sys.version}),flush=True)
 with ProcessPoolExecutor(max_workers=min(15,a.workers)) as ex:
  fs={ex.submit(task,j,barrier if a.workers>=15 else None):j for j in range(15)}
  for f in as_completed(fs):
   j=fs[f]
   try:
    r=f.result();done.append(r);failed+=r['bad'];print(json.dumps({'event':'done','job':j,'pid':r['worker_pid'],'cache':r['cache_state'],'evals':r['evaluations'],'valid':r['run_validation']}),flush=True)
   except Exception as e:failed.append({'job':j,'error':repr(e)});print(json.dumps({'event':'fail','job':j,'error':repr(e)}),flush=True)
 manager.shutdown();best=min((r['best'] for r in done),key=lambda z:z['value']);laws=set(s for r in done for s in r['law_signatures']);out={'version':VERSION,'evaluations':sum(r['evaluations'] for r in done),'unique_probability_laws':len(laws),'failures':failed,'best':best,'records':sorted(done,key=lambda z:z['job']),'seconds':time.time()-start,'full_pool':a.workers>=15}
 (ROOT/'search_failed_inputs.json').write_text(json.dumps(failed,indent=2));(ROOT/'search_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps({'event':'finish','evaluations':out['evaluations'],'failed':len(failed),'seconds':out['seconds']}),flush=True)
 if failed:raise SystemExit(1)
if __name__=='__main__':main()
