"""Finite floating-point checks of the strict face-mass fibres.
General all-set and Brunn-Minkowski arguments are analytic, not certified by this code.
"""
import concurrent.futures
import itertools
import json
import math
import os
import platform
import scipy
from pathlib import Path
import time
import warnings
import numpy as np
from scipy.optimize import linprog,OptimizeWarning
from scipy.spatial import HalfspaceIntersection,ConvexHull


def inequalities(n,beta,f,budget=1.):
    e=sum(f);rows=[]
    for mask in itertools.product([0,1],repeat=2*n):
        q=n-beta-sum(mask)/2
        u=-np.array([mask[2*j]-mask[2*j+1] for j in range(n)],float)
        rhs=budget+(1-beta/n)*e-np.dot(mask,f)
        rows.append((q,u,rhs))
    positive=[r for r in rows if r[0]>1e-12]
    negative=[r for r in rows if r[0]<-1e-12]
    output=[(u,b) for q,u,b in rows if abs(q)<=1e-12]
    for qp,up,bp in positive:
        for qn,un,bn in negative:
            output.append((qp*un-qn*up,qp*bn-qn*bp))
    unique={}
    for a,b in output:
        norm=np.linalg.norm(a)
        if norm<1e-12:
            # A zero normal with zero RHS forbids STRICT feasibility.
            if b<=1e-10:
                return None,None,'strict-zero-normal-obstruction'
            continue
        a=a/norm;b=b/norm;key=tuple(np.round(a,13))
        if key not in unique or b<unique[key][1]:
            unique[key]=(a,b)
    return np.array([v[0] for v in unique.values()]),np.array([v[1] for v in unique.values()]),''


def fibre_volume(n,beta,f):
    a,b,note=inequalities(n,beta,np.asarray(f))
    if a is None:
        return 0.,note
    # HiGHS receives threads=1; SciPy warns that this solver-specific option is passed through.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',category=OptimizeWarning,message='Unrecognized options detected.*')
        lp=linprog(np.r_[np.zeros(n),-1],A_ub=np.column_stack([a,np.ones(len(a))]),
                   b_ub=b,bounds=[(None,None)]*(n+1),method='highs-ds',options={'threads':1})
    assert lp.success,lp.message
    if lp.x[-1]<=1e-9:
        return 0.,'no-positive-projected-interior-margin'
    vertices=HalfspaceIntersection(np.column_stack([a,-b]),lp.x[:n]).intersections
    volume=float(ConvexHull(vertices).volume)
    assert np.max(a@vertices.T-b[:,None])<2e-8
    return volume,'positive-strict-margin'


def leading(n,beta):
    c=min(beta,n-beta);k=math.floor(2*c)
    return 2**n/(math.factorial(k)*(2*c)**(n-k))


def run_parameter(job):
    n,beta=job;rng=np.random.default_rng(20261105+100*n+int(100*beta))
    shape=np.linspace(-1,1,n)
    hs=[np.zeros(n),*[a*shape for a in [.1,.25,.49,.5,1.,2.]]]
    hs += [h-h.mean() for h in rng.uniform(-1.25,1.25,size=(20,n))]
    top=leading(n,beta);ratios=[];invariance=[];zero=0;worst=None
    for h in hs:
        f=np.repeat(h,2);volume,note=fibre_volume(n,beta,f)
        zero+=int(volume==0);ratios.append(volume/top)
        if worst is None or volume/top>worst['ratio']:
            worst={'f':f.tolist(),'ratio':volume/top}
        shift=rng.uniform(-.7,.7,n);scale=float(rng.uniform(-.7,.7))
        moved=f+scale/2+np.ravel(np.column_stack([shift,-shift]))
        other,_=fibre_volume(n,beta,moved)
        invariance.append(abs(volume-other)/max(1.,top))
    assert max(ratios)<1+2e-8
    assert max(invariance)<2e-8
    return {'n':n,'beta':beta,'base_cases':len(hs),'volume_evaluations':2*len(hs),
            'max_volume_over_L':max(ratios),'empty_base_cases':zero,
            'max_translation_scale_invariance_error':max(invariance),'worst_case':worst}


def cross_cases():
    out=[]
    for n,beta in [(2,.5),(2,1.),(2,1.5),(3,.5),(3,1.5),(3,2.5)]:
        for amplitude in [0.,.1,.25,.5]:
            h=amplitude*np.linspace(-1,1,n)
            volume,note=fibre_volume(n,beta,np.repeat(h,2))
            out.append({'n':n,'beta':beta,'amplitude':amplitude,'volume':volume,'status':note})
    return out


def compare(base,rows):
    matlab=base/'results/global_leading_matlab.json'
    if not matlab.exists():
        return
    differences=[]
    for n,beta,amplitude,volume,*_ in json.loads(matlab.read_text())['rows']:
        ref=next(r for r in rows if r['n']==n and r['beta']==beta and r['amplitude']==amplitude)
        differences.append(abs(volume-ref['volume'])/max(1.,leading(n,beta)))
    result={'kind':'cross-tool-numerical-comparison','cases':len(differences),
            'max_normalized_volume_difference':max(differences)}
    (base/'results/global_leading_crosscheck.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


def main():
    import sys
    base=Path(__file__).resolve().parents[1]
    if '--compare-only' in sys.argv:
        compare(base,json.loads((base/'results/global_leading_python.json').read_text())['cross_cases'])
        return
    started=time.perf_counter()
    params=[(2,b) for b in [.25,.5,.75,1.,1.25,1.5,1.75]]+[(3,b) for b in [.5,1.,1.5,2.,2.5]]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        summary=list(pool.map(run_parameter,params))
    critical=[]
    for delta in [0.,.25,.49,.5-1e-6,.5,.5+1e-6,.75]:
        volume,note=fibre_volume(2,1.,[delta,delta,-delta,-delta])
        expected=2. if 2*delta<1 else 0.
        assert abs(volume-expected)<2e-9
        critical.append({'delta':delta,'strict_volume':volume,'expected':expected,'status':note})
    cross=cross_cases()
    result={'kind':'non-certified-finite-polytope-numerical-evidence','seed_base':20261105,
            'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
            'precision':'IEEE754 double','workers':4,'logical_cpus':os.cpu_count(),
            'HiGHS_threads_per_worker':1,'elapsed_seconds':time.perf_counter()-started,
            'zero_normal_tolerance':1e-12,'strict_RHS_guard':1e-10,'interior_radius_guard':1e-9,
            'coverage':'Finite deterministic and seeded face vectors; no general optimization certificate.',
            'parameter_summaries':summary,'critical_strict_threshold':critical,'cross_cases':cross}
    (base/'results/global_leading_python.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Volume computations:',sum(r['volume_evaluations'] for r in summary)+len(cross)+len(critical))
    print('Max volume/L:',max(r['max_volume_over_L'] for r in summary))
    print('Max normalized symmetry residual:',max(r['max_translation_scale_invariance_error'] for r in summary))
    print('Elapsed:',result['elapsed_seconds'])
    compare(base,cross)

if __name__=='__main__':
    main()
