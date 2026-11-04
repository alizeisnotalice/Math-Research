#!/usr/bin/env python3
"""Minimal 120-digit audit of bounded truncations of the cube angular rearrangement.

The continuous angular law has density (n-1)/(n gamma)*a^(-1/gamma-1).
The exponent in its q-moment is q-1/gamma (not gamma*q-1).
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import argparse, hashlib, json, os as _os, platform, shutil, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from decimal import Decimal as D, getcontext
from fractions import Fraction
from multiprocessing import Manager
from pathlib import Path

ROOT=Path(__file__).resolve().parent; CACHE=ROOT/'runtime_cache'; SUMMARY=ROOT/'round29_summary.json'
VERSION='round29-archive-v2-20261105'; PREC=120
NS=(2,3); GS=('0.25','0.5','0.75'); CS=('0.25','0.5','0.75','0.95'); RS=(2,4,16,256)
def source_hash(): return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def pw(x,e): return D(0) if x==0 else (e*x.ln()).exp()
def full(n,g,q): return (D(n)-g*q)/(D(n)*(1-g*q))
def partial(n,g,q,A):
    w=D(n-1)/n; e=q-D(1)/g; z=g*q-1
    return D(1)/n+(w/g*A.ln() if e==0 else w*(pw(A,e)-1)/z)
def tail(n,g,q,A): return D(n-1)/n*pw(A,q-D(1)/g)/(1-g*q)
def exact_rational_checks():
    # Both choices make q-1/gamma an integer, so this is exact Fraction arithmetic.
    tests=((Fraction(0),Fraction(2)),(Fraction(1),Fraction(4)))
    for q,A in tests:
        n=2; g=Fraction(1,2); w=Fraction(n-1,n); e=q-1/g
        le=Fraction(1,n)+w*(A**e-1)/(g*q-1)
        ta=w*A**e/(1-g*q); fu=(Fraction(n)-g*q)/(n*(1-g*q))
        assert le+ta==fu
    getcontext().prec=PREC; gdec=D('0.5')
    # At q=1/gamma only the finite partial moment is used: 1/2+log(2).
    assert abs(partial(2,gdec,D(2),D(2))-(D(1)/2+D(2).ln())) < D('1e-100')
    return ['n2_g1_2_q0_A2','n2_g1_2_q1_A4','n2_g1_2_q2_A2_log_partial']
def record(item):
    n,gs,cs=item; getcontext().prec=PREC; g=D(gs); c=D(cs); th=(1-g)*c; t=D(n)*th; p=D(1)/(1-th)
    ea=full(n,g,D(1)); ep=full(n,g,p); ell=pw(ep,D(1)/p)
    for q in (D(0),D(1),p):
        for A in (D(1),D(2),D(257)):
            assert abs(partial(n,g,q,A)+tail(n,g,q,A)-full(n,g,q)) < D('1e-100')
    vals=[]
    for R0 in RS:
        R=D(R0); A=ell*pw(R,D(n)-t); lep=partial(n,g,p,A)
        v=lep/ep+pw(R,n)*tail(n,g,D(0),A)
        m=pw(ell,-t/(D(n)-t))*lep+pw(R,t)*tail(n,g,D(1),A)
        rel=m*pw(v,-th)/ell; assert D(0)<v<=1 and D(0)<rel<=1
        vals.append((R0,rel))
    return {'n':n,'gamma':gs,'c':cs,'R256_relative':str(vals[-1][1]),'max_deficit':str(max(1-x for _,x in vals)),'pid':_os.getpid()}
def signature(item):
    p={'version':VERSION,'source_sha256':source_hash(),'item':item,'precision':PREC,'python':sys.version,'platform':platform.platform()}
    return hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest()
def job(item,barrier=None):
    name='n%d_g%s_c%s.json'%(item[0],item[1].replace('.','p'),item[2].replace('.','p')); path=CACHE/name; sig=signature(item)
    if barrier is not None: barrier.wait(timeout=60)
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('signature')==sig: return old|{'cache':'reused','replay_pid':_os.getpid()}
    out=record(item)|{'signature':sig,'cache':'computed'}; tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out));tmp.replace(path);return out
def compact(rows):
    wanted={(2,'0.5','0.5'),(2,'0.75','0.95'),(3,'0.5','0.5'),(3,'0.75','0.95')}
    return [{k:v for k,v in r.items() if k in ('n','gamma','c','R256_relative','max_deficit')} for r in rows if (r['n'],r['gamma'],r['c']) in wanted]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=('fresh','replay'));a=ap.parse_args(); items=[(n,g,c) for n in NS for g in GS for c in CS]
    rational_checks=exact_rational_checks() # Always execute; do not let an old summary bypass it.
    if a.mode=='fresh' and CACHE.exists(): shutil.rmtree(CACHE)
    CACHE.mkdir(exist_ok=True); rows=[]; bad=[]
    with Manager() as mgr:
      barrier=mgr.Barrier(15)
      with ProcessPoolExecutor(max_workers=15) as ex:
        fs={ex.submit(job,x,barrier):x for x in items[:15]}
        for f in as_completed(fs):
          try: rows.append(f.result())
          except Exception as e: bad.append({'item':fs[f],'error':repr(e)})
        fs={ex.submit(job,x):x for x in items[15:]}
        for f in as_completed(fs):
          try: rows.append(f.result())
          except Exception as e: bad.append({'item':fs[f],'error':repr(e)})
    rows.sort(key=lambda r:(r['n'],D(r['gamma']),D(r['c']))); assert not bad
    data=json.loads(SUMMARY.read_text()) if SUMMARY.exists() else {'version':VERSION,'grid':{'n':[2,3],'gamma':['.25','.5','.75'],'c':['.25','.5','.75','.95'],'R':list(RS)}}
    data['rational_partition_checks']=rational_checks
    data[a.mode]={'tasks':len(rows),'failures':len(bad),'worker_pids':sorted(set(r['pid'] for r in rows)),'computed':sum(r['cache']=='computed' for r in rows),'reused':sum(r['cache']=='reused' for r in rows),'samples':compact(rows)}
    SUMMARY.write_text(json.dumps(data,indent=2)); print(json.dumps(data[a.mode],indent=2))
if __name__=='__main__': main()
