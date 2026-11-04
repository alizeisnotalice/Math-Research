#!/usr/bin/env python3
"""15-process candidate search for binary mean-one L1 coefficients.

For a in (0,1), probability q of a, and b=(1-qa)/(1-q), evaluate
c_p=(E f^p-1)/(E|f-1|)^2 using the centered Bregman numerator.  All output
is numerical candidate evidence only; p=3 has a separate exact audit.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[_k]='1'
import argparse, hashlib, json, math, os as _os, platform, sys, time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path

ROOT=Path(__file__).resolve().parent; CHECK=ROOT/'checkpoints'
VERSION='binary-l1-v1-bregman-grid-20261105'
PS=(1.1,1.5,2.0,3.0,4.0,8.0)

def bregman(x,p): return math.expm1(p*math.log(x))-p*(x-1.0)
def ratio(q,a,p):
    delta=q*(1-a); b=(1-q*a)/(1-q)
    u=q*bregman(a,p)+(1-q)*bregman(b,p)
    return u/(4*delta*delta),b,u

def grid_for_p(p):
    # Boundary probes, an ordinary coarse grid, and p=3 local refinement.
    qs={i/200 for i in range(1,200)}
    qs.update(2.0**(-k) for k in range(2,21));qs.update(1-2.0**(-k) for k in range(2,21))
    if p==3.0:
        qstar=(3-math.sqrt(5))/2
        qs.update(qstar+j*1e-5 for j in range(-300,301) if 0<qstar+j*1e-5<1)
    aa={i/100 for i in range(1,100)}
    aa.update(2.0**(-k) for k in range(1,41))
    return sorted(qs),sorted(aa)

def signature(job):
    return hashlib.sha256(json.dumps({'job':job,'p':PS,'grid':'q200+dyadic+P3refine;a100+dyadic40'},sort_keys=True).encode()).hexdigest()

def task(job,barrier=None):
    if barrier is not None: barrier.wait(timeout=30)
    path=CHECK/f'job_{job:02d}.json';sig=signature(job)
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('version')==VERSION and old.get('signature')==sig:
            old['cache_state']='reused';old['checkpoint_worker_pid']=old.get('worker_pid');old['worker_pid']=_os.getpid();return old
    start=time.time(); best={}; count=0; bad=[]
    for p in PS:
        qs,aa=grid_for_p(p)
        for qi,q in enumerate(qs):
            for ai,a in enumerate(aa):
                if (qi*len(aa)+ai)%15 !=job: continue
                try:
                    value,b,u=ratio(q,a,p); count+=1
                    if not (math.isfinite(value) and value>0 and b>1 and u>0): bad.append({'p':p,'q':q,'a':a,'value':value})
                    elif p not in best or value<best[p]['value']: best[p]={'value':value,'q':q,'a':a,'b':b,'bregman_u':u}
                except Exception as exc: bad.append({'p':p,'q':q,'a':a,'exception':repr(exc)})
    out={'job':job,'version':VERSION,'signature':sig,'worker_pid':_os.getpid(),'checkpoint_worker_pid':_os.getpid(),
         'cache_state':'computed','seconds':time.time()-start,'evaluations':count,'best':best,'bad':bad,'run_validation':not bad}
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2,allow_nan=False));tmp.replace(path);return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);args=ap.parse_args();CHECK.mkdir(exist_ok=True)
    from multiprocessing import Manager
    manager=Manager(); barrier=manager.Barrier(15); start=time.time(); done=[]; failed=[]
    print(json.dumps({'event':'start','tasks':15,'workers':args.workers,'version':VERSION,'platform':platform.platform(),'python':sys.version}),flush=True)
    with ProcessPoolExecutor(max_workers=min(15,args.workers)) as ex:
        fs={ex.submit(task,j,barrier if args.workers>=15 else None):j for j in range(15)}
        for f in as_completed(fs):
            j=fs[f]
            try:
                r=f.result();done.append(r)
                if not r['run_validation']: failed += r['bad']
                print(json.dumps({'event':'done','job':j,'pid':r['worker_pid'],'cache':r['cache_state'],'evals':r['evaluations'],'valid':r['run_validation']}),flush=True)
            except Exception as exc:
                failed.append({'job':j,'exception':repr(exc)});print(json.dumps({'event':'fail','job':j,'error':repr(exc)}),flush=True)
    manager.shutdown(); merged={}
    for r in done:
        for pk,z in r['best'].items():
            p=float(pk)
            if p not in merged or z['value']<merged[p]['value']: merged[p]=z
    summary={'version':VERSION,'tasks':len(done),'evaluations':sum(r['evaluations'] for r in done),'failures':failed,
             'best':merged,'records':sorted(done,key=lambda x:x['job']),'seconds':time.time()-start,'full_pool':args.workers>=15}
    (ROOT/'failed_inputs.json').write_text(json.dumps(failed,indent=2,allow_nan=False))
    (ROOT/'binary_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    print(json.dumps({'event':'finish','tasks':len(done),'evaluations':summary['evaluations'],'failed':len(failed),'seconds':summary['seconds']}),flush=True)
    if failed: raise SystemExit(1)
if __name__=='__main__':main()
