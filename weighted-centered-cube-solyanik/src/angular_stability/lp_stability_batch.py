#!/usr/bin/env python3
"""Finite-probability diagnostics for the stated Lp stability inequalities.

This is deliberately an angular/probability-vector test, not an outer-halo or
maximal-average computation.  Inputs are rational positive weights and values
with exact weighted mean one.  Float checks are diagnostics only; integer-p
rows additionally carry exact Fraction checks.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse, hashlib, json, math, os as _os, platform, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parent; CHECK=ROOT/'checkpoints'
VERSION='lp-angular-v2-sharp-spike-scale10-20261105'
PS=(Q(11,10),Q(5,4),Q(3,2),Q(19,10),Q(2),Q(21,10),Q(3),Q(4),Q(8),Q(16))
INTEGER_PS={2,3,4,8,16}

def norm_mean(weights, values):
    mean=sum((r*x for r,x in zip(weights,values)),Q(0))
    return weights,[x/mean for x in values]

def vec_balanced(delta):
    return [Q(1,2),Q(1,2)],[1-delta,1+delta]

def vec_hole(q):
    return [q,1-q],[q,1+q]

def vec_spike(p,scale):
    # M=2^(den(p)*scale), eps=M^(-(p+1)) is rational exactly.
    M=Q(2**(p.denominator*scale),1)
    exponent=(p.numerator+p.denominator)*scale
    eps=Q(1,2**exponent)
    c=(1-eps*M)/(1-eps)
    assert eps>0 and c>0 and eps*M+(1-eps)*c==1
    return [eps,1-eps],[M,c]

def vec_random(seed):
    # Deterministic positive rational angular data; common denominator is exact.
    raw_w=[Q((seed*7+3*i)%11+1,37) for i in range(4)]
    weights=[x/sum(raw_w,Q(0)) for x in raw_w]
    raw_x=[Q((seed*13+5*i)%29+1,(seed+i)%7+2) for i in range(4)]
    return norm_mean(weights,raw_x)

def input_vectors(seed,p):
    out=[]
    for power in (3,8,16): out.append((f'balanced_2^-{power}',*vec_balanced(Q(1,2**power))))
    for power in (3,6,9): out.append((f'hole_2^-{power}',*vec_hole(Q(1,2**power))))
    for scale in (1,3,10): out.append((f'spike_scale_{scale}',*vec_spike(p,scale)))
    out.append(('random_angular',*vec_random(seed)))
    return out

def bregman_float(f,p):
    d=f-1.0
    return math.expm1(p*math.log1p(d))-p*d

def float_metrics(weights,values,p):
    rf=[float(r) for r in weights]; xf=[float(x) for x in values]
    # Since Ef=1, this is E[f^p-1-p(f-1)]=E f^p-1 with no cancellation.
    u=math.fsum(r*bregman_float(x,p) for r,x in zip(rf,xf))
    dp=math.fsum(r*abs(x-1.0)**p for r,x in zip(rf,xf))
    d2=math.fsum(r*(x-1.0)**2 for r,x in zip(rf,xf))
    if p<=2:
        rhs=math.expm1((2.0/p)*math.log1p(u))/(p-1.0)
        lhs=dp**(2.0/p)
        return {'bregman_u':u,'lhs':lhs,'rhs':rhs,'slack':rhs-lhs,
                'sharp_ratio':lhs/rhs if rhs else float('nan'),'roundoff_note':'Bregman/expm1/log1p'}
    rhs2=u/(p-1.0)
    return {'bregman_u':u,'lp_power':dp,'lp_rhs':u,'lp_slack':u-dp,
            'variance':d2,'variance_rhs':rhs2,'variance_slack':rhs2-d2,
            'spike_ratio':dp/u if u else float('nan'),
            'hole_ratio':(p-1.0)*d2/u if u else float('nan'),'roundoff_note':'Bregman centered'}

def exact_integer_metrics(weights,values,p):
    n=int(p); mean=sum((r*x for r,x in zip(weights,values)),Q(0)); assert mean==1
    b=sum((r*(x**n-1-n*(x-1)) for r,x in zip(weights,values)),Q(0))
    dp=sum((r*abs(x-1)**n for r,x in zip(weights,values)),Q(0))
    d2=sum((r*(x-1)**2 for r,x in zip(weights,values)),Q(0))
    scalar=all(abs(x-1)**n <= x**n-1-n*(x-1) and (n-1)*(x-1)**2 <= x**n-1-n*(x-1) for x in values)
    if n==2:
        valid=(dp==b and d2==b and scalar)
    else:
        valid=(dp<=b and (n-1)*d2<=b and scalar)
    return {'exact_bregman':str(b),'exact_lp_power':str(dp),'exact_variance':str(d2),'exact_scalar_pass':scalar,'exact_averaged_pass':valid}

def row(seed,p,label,w,x):
    pf=float(p); fm=float_metrics(w,x,pf)
    # Tolerance only classifies arithmetic diagnostics, never proof.
    ok=(fm['slack']>=-2e-12*max(1.0,abs(fm['lhs']),abs(fm['rhs'])) if pf<=2 else
        fm['lp_slack']>=-2e-12*max(1.0,abs(fm['lp_rhs'])) and fm['variance_slack']>=-2e-12*max(1.0,abs(fm['variance_rhs'])))
    result={'seed':seed,'p':str(p),'label':label,'weights':[str(z) for z in w],'values':[str(z) for z in x],
            'float':fm,'float_diagnostic_pass':ok,'exact':None}
    if p.denominator==1 and int(p) in INTEGER_PS:
        result['exact']=exact_integer_metrics(w,x,p)
        result['run_validation']=ok and result['exact']['exact_averaged_pass']
    else: result['run_validation']=ok
    return result

def signature(seed):
    payload={'seed':seed,'p':[str(p) for p in PS],'families':'balanced3,hole3,spike1,3,10,random1'}
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()

def task(seed,barrier=None):
    if barrier is not None: barrier.wait(timeout=30)
    path=CHECK/f'group_{seed:02d}.json'; sig=signature(seed)
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('version')==VERSION and old.get('signature')==sig:
            old['cache_state']='reused'; old['checkpoint_worker_pid']=old.get('worker_pid'); old['worker_pid']=_os.getpid(); return old
    start=time.time(); rows=[]
    for p in PS:
        for label,w,x in input_vectors(seed,p): rows.append(row(seed,p,label,w,x))
    failed=[r for r in rows if not r['run_validation']]
    output={'seed':seed,'version':VERSION,'signature':sig,'worker_pid':_os.getpid(),'checkpoint_worker_pid':_os.getpid(),
            'cache_state':'computed','seconds':time.time()-start,'rows':rows,'run_validation':not failed,'failed':failed}
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(output,indent=2,allow_nan=False));temp.replace(path);return output

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);args=ap.parse_args();CHECK.mkdir(exist_ok=True)
    from multiprocessing import Manager
    manager=Manager();barrier=manager.Barrier(15);out=[];failure=[];start=time.time()
    print(json.dumps({'event':'start','workers':args.workers,'tasks':15,'version':VERSION,'platform':platform.platform(),'python':sys.version}),flush=True)
    with ProcessPoolExecutor(max_workers=min(15,args.workers)) as ex:
        fs={ex.submit(task,s,barrier if args.workers>=15 else None):s for s in range(15)}
        for future in as_completed(fs):
            seed=fs[future]
            try:
                r=future.result();out.append(r)
                if not r['run_validation']: failure.extend(r['failed'])
                print(json.dumps({'event':'done','seed':seed,'pid':r['worker_pid'],'cache':r['cache_state'],'valid':r['run_validation']}),flush=True)
            except Exception as exc:
                failure.append({'seed':seed,'exception':repr(exc)})
                print(json.dumps({'event':'fail','seed':seed,'error':repr(exc)}),flush=True)
    manager.shutdown()
    summary={'version':VERSION,'tasks':len(out),'failures':len(failure),'seconds':time.time()-start,
             'records':sorted(out,key=lambda z:z['seed']),'full_pool':args.workers>=15}
    (ROOT/'failed_inputs.json').write_text(json.dumps(failure,indent=2,allow_nan=False))
    (ROOT/'lp_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    print(json.dumps({'event':'finish','tasks':len(out),'failed':len(failure),'seconds':summary['seconds']}),flush=True)
    if failure: raise SystemExit(1)

if __name__=='__main__': main()
