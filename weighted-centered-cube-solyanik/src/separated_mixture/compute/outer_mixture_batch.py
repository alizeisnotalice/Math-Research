#!/usr/bin/env python3
"""Outer-weight mixture experiments: Lebesgue halo, weighted cost afterwards.

This is intentionally separate from rejected/mixture_halo_batch.py.  The halo
H_alpha(E) is computed only from |E intersection [x-r,x+r]|/(2r)>alpha.
Weights enter solely through w(H)/w(E).  `halo_fraction` is a local extraction
of experiments/round013/verify.py; `halo_float_cells` is independently adapted
from experiments/round013/search.py's float piecewise-linear implementation.
Both exploit the same one-dimensional piecewise-linear geometry, so agreement is
a regression check rather than two independent mathematical proofs.
"""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import argparse, json, math, platform, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parent
CHECKPOINTS=ROOT/'outer_checkpoints'

def merge(intervals):
    out=[]
    for a,b in sorted((a,b) for a,b in intervals if b>a):
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out

def halo_fraction(E,alpha):
    # Reused, without import side effects, from round013/verify.py.
    ends=sum(([a,b] for a,b in E),[]); ell=sum(b-a for a,b in E)
    lo=ends[0]-ell/(2*alpha); hi=ends[-1]+ell/(2*alpha); pieces=list(E)
    for e in ends:
        breaks=sorted({lo,hi,*((e+f)/2 for f in ends if lo<(e+f)/2<hi)})
        def G(x):
            left,right=min(e,2*x-e),max(e,2*x-e)
            return sum(max(F(0),min(right,b)-max(left,a)) for a,b in E)-2*alpha*abs(x-e)
        for a,b in zip(breaks,breaks[1:]):
            ga,gb=G(a),G(b)
            if ga>0 and gb>0: pieces.append((a,b))
            elif ga>0: pieces.append((a,a+(b-a)*ga/(ga-gb)))
            elif gb>0: pieces.append((a+(b-a)*(-ga)/(gb-ga),b))
    return merge(pieces)

def halo_float_cells(E,alpha):
    # Separate floating endpoint-cell implementation based on round013/search.py.
    E=[(float(a),float(b)) for a,b in E]; alpha=float(alpha)
    ends=[v for ab in E for v in ab]; ell=sum(b-a for a,b in E)
    lo=ends[0]-ell/(2*alpha); hi=ends[-1]+ell/(2*alpha); pieces=list(E)
    def mass(a,b): return sum(max(0.,min(b,d)-max(a,c)) for c,d in E)
    for e in ends:
        cuts=sorted({lo,hi,*((e+v)/2 for v in ends if lo<(e+v)/2<hi)})
        for a,b in zip(cuts,cuts[1:]):
            def g(x): return mass(min(e,2*x-e),max(e,2*x-e))-2*alpha*abs(x-e)
            ga,gb=g(a),g(b); eps=1e-13
            if abs(ga)<eps: ga=0.
            if abs(gb)<eps: gb=0.
            if ga>0 and gb>0: pieces.append((a,b))
            elif ga>0: pieces.append((a,a+(b-a)*ga/(ga-gb)))
            elif gb>0: pieces.append((a+(b-a)*(-ga)/(gb-ga),b))
    return merge(pieces)

def signed_power(x,tau): return math.copysign(abs(x)**tau,x) if x else 0.

def component_mass(component,a,b):
    if component[0]=='base': return component[1]*(b-a)
    _,c,z,tau=component
    return c*(signed_power(b-z,tau)-signed_power(a-z,tau))/tau

def mass(component,intervals): return sum(component_mass(component,float(a),float(b)) for a,b in intervals)

def configs():
    taus=(.15,.25,.5,.75,.9); result=[]
    Es=[[(F(-6,5),F(-1,20)),(F(1,20),F(6,5))],
        [(F(-1),F(-1,10)),(F(1,15),F(9,10))],
        [(F(-7,8),F(-1,8)),(F(1,10),F(19,20))]]
    for k in range(15):
        if k>=11: z1,z2=F(-3-k%3,25),F(3+k%3,25)
        else: z1,z2=F(-18+2*(k%4),25),F(10+2*((3*k)%4),25)
        result.append({'id':k,'c0':.2+.05*(k%3),
          'singular':[(.35+.13*(k%5),z1,taus[k%5]),(.55+.11*((2*k)%5),z2,taus[(2*k+1)%5])],
          'E':Es[k%3]})
    return result

def local_rows(config,j):
    c,z,tau=config['singular'][j]; other=config['singular'][1-j][1]; sep=abs(z-other)
    comps=[('base',config['c0'])]+[('singular',float(cc),float(zz),tt) for cc,zz,tt in config['singular']]
    alpha=F(9,10); rows=[]
    for k in range(4,11):
        eps=F(1,2**k)
        if 5*eps>=sep: continue
        # A two-interval input with a thin hole containing z; alpha=.9 sees the hole.
        E=[(z-F(3,2)*eps,z-F(1,20)*eps),(z+F(1,20)*eps,z+F(3,2)*eps)]
        H=halo_fraction(E,alpha)
        den=mass(comps[j+1],E); pure=mass(comps[j+1],H)/den
        mixed_den=sum(mass(q,E) for q in comps); mixed=sum(mass(q,H) for q in comps)/mixed_den
        scale=float(eps)**(1-tau)
        rows.append({'k':k,'epsilon':float(eps),'pure':pure,'mixed':mixed,
                     'difference':mixed-pure,'scaled_abs_difference':abs(mixed-pure)/scale,
                     'halo_intervals':len(H)})
    return {'component':j,'site':str(z),'tau':tau,'separation':str(sep),'rows':rows}

def task(config, start_barrier=None):
    # The barrier is an audit device: with exactly 15 submitted tasks it proves
    # that the stated 15-process run actually reached 15 live workers.
    if start_barrier is not None:
        start_barrier.wait(timeout=30)
    file=CHECKPOINTS/f"config_{config['id']:02d}.json"
    expected_config={'id':config['id'],'c0':config['c0'],
                     'singular':[[c,str(z),t] for c,z,t in config['singular']],
                     'E':[[str(a),str(b)] for a,b in config['E']]}
    if file.exists():
        cached=json.loads(file.read_text())
        if (cached.get('parameter_version')=='outer-v2-20261105'
            and cached.get('config')==expected_config
            and len(cached.get('rows',[]))==5
            and all(r['float_agrees'] and r['convex_pass'] for r in cached['rows'])):
            cached['reused']=True
            cached['current_worker_pid']=os.getpid()
            return cached
    start=time.time(); E=config['E']
    comps=[('base',config['c0'])]+[('singular',float(c),float(z),tau) for c,z,tau in config['singular']]
    alphas=[F(1,4),F(2,5),F(3,4),F(9,10),F(99,100)]
    rows=[]
    for alpha in alphas:
        H=halo_fraction(E,alpha); HF=halo_float_cells(E,alpha)
        if len(H)==len(HF):
            endpoint_differences=[]
            for (left,right),(left_float,right_float) in zip(H,HF):
                endpoint_differences.extend((abs(float(left)-left_float),abs(float(right)-right_float)))
            error=max(endpoint_differences,default=0.)
        else:
            error=float('inf')
        ratios=[mass(q,H)/mass(q,E) for q in comps]
        weights=[mass(q,E) for q in comps]
        mixture=sum(mass(q,H) for q in comps)/sum(weights)
        rows.append({'alpha':str(alpha),'halo_intervals':len(H),'float_agrees':error<2e-10,
                     'endpoint_error':error,'component_common_halo_ratios':ratios,
                     'mixture_ratio':mixture,'convex_upper':max(ratios),
                     'convex_pass':mixture<=max(ratios)+2e-12})
    mapping=[]
    for m in (2,3,4,8):
        root=F(m-1,m)
        for side, rho in [('below',root-F(1,1000)),('above',root+F(1,1000))]:
            lam=rho**m; phi=1-m*(1-rho)
            mapping.append({'m':m,'side':side,'lambda':str(lam),'phi':str(phi),
                            'halo_defined':phi>0})
    rec={'config':{'id':config['id'],'c0':config['c0'],'singular':[[c,str(z),t] for c,z,t in config['singular']],
                   'E':[[str(a),str(b)] for a,b in E]},'rows':rows,'critical_mapping':mapping,
         'local_blowup':[local_rows(config,j) for j in (0,1)],'seconds':time.time()-start,'status':'completed',
         'worker_pid':os.getpid(),'current_worker_pid':os.getpid(),'reused':False,'python_version':sys.version,'parameter_version':'outer-v2-20261105'}
    tmp=file.with_suffix('.tmp'); tmp.write_text(json.dumps(rec,indent=2,allow_nan=False));tmp.replace(file); return rec

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);ap.add_argument('--audit-full-pool',action='store_true');a=ap.parse_args();CHECKPOINTS.mkdir(exist_ok=True)
    start=time.time();cs=configs();done=[];fails=[]
    print(json.dumps({'event':'start','tasks':15,'workers':min(a.workers,15),'python':sys.version,'platform':platform.platform()}),flush=True)
    barrier=None; manager=None
    if a.audit_full_pool and len(cs)==15 and min(a.workers,15)==15:
        from multiprocessing import Manager
        manager=Manager();barrier=manager.Barrier(15)
    with ProcessPoolExecutor(max_workers=min(a.workers,15)) as ex:
      fs={ex.submit(task,c,barrier):c['id'] for c in cs}
      for f in as_completed(fs):
        try:
          r=f.result();done.append(r);print(json.dumps({'event':'done','id':fs[f],'pid':r['worker_pid'],'seconds':r['seconds'],
             'outer_geometry_all_agree':all(x['float_agrees'] for x in r['rows']),
             'convex_all_pass':all(x['convex_pass'] for x in r['rows'])}),flush=True)
        except Exception as e: fails.append({'id':fs[f],'error':repr(e)});print(json.dumps({'event':'fail','id':fs[f],'error':repr(e)}),flush=True)
    out={'status':'completed' if not fails else 'completed_with_failures','finished':len(done),'failures':fails,
         'seconds':time.time()-start,'audit_full_pool':bool(a.audit_full_pool),'records':sorted(done,key=lambda x:x['config']['id'])}
    (ROOT/'outer_mixture_summary.json').write_text(json.dumps(out,indent=2,allow_nan=False));print(json.dumps({'event':'finish','finished':len(done),'failures':len(fails),'seconds':out['seconds']}),flush=True)
    if manager is not None: manager.shutdown()
    valid=(not fails and len(done)==15
           and all(r['float_agrees'] and r['convex_pass'] for rec in done for r in rec['rows'])
           and all(math.isfinite(r['mixed']) and math.isfinite(r['pure'])
                   for rec in done for part in rec['local_blowup'] for r in part['rows']))
    print(json.dumps({'validation':'passed' if valid else 'failed',
                      'computed_tasks':sum(not r.get('reused',False) for r in done),
                      'reused_tasks':sum(r.get('reused',False) for r in done),
                      'current_worker_pids':sorted({r['current_worker_pid'] for r in done})}),flush=True)
    if not valid: raise SystemExit(1)
if __name__=='__main__': main()
