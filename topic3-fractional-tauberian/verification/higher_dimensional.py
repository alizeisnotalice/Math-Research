"""Deterministic Gauss quadrature of all contact roots; numerical evidence only.
Root implementation uses explicit sorted-contact products; MATLAB independently
uses the original minimum-of-three interval overlap and operator bisection.
"""
import concurrent.futures
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import resource
import time
import sys
import numpy as np
from scipy.special import roots_legendre
import scipy


def directions(n,q):
    if n==1:
        return np.ones((1,1)),np.ones(1)
    nodes,weights=roots_legendre(q)
    nodes=(nodes+1)/2;weights=weights/2
    mesh=np.meshgrid(*([nodes]*(n-1)),indexing='ij')
    wm=np.meshgrid(*([weights]*(n-1)),indexing='ij')
    u=np.column_stack([x.ravel() for x in mesh])
    w=np.prod(np.column_stack([x.ravel() for x in wm]),axis=1)
    w*=np.prod(u**np.arange(n-2,-1,-1)[None,:],axis=1)
    z=np.column_stack([np.ones(len(u)),np.cumprod(u,axis=1)])
    assert abs(w.sum()-1/math.factorial(n-1))<2e-14
    return z,w


def contact_roots(z,beta,t):
    n=z.shape[1];sum_prev=np.cumsum(z,axis=1)-z;idx=np.arange(n)
    am=sum_prev+(2*beta-idx)*z
    ap=sum_prev+(2*(n-beta)-idx)*z
    answer=[];worst=0.
    for sign in [-1,1]:
        for i in range(n):
            zi=z[:,i];r=t/(am[:,i] if sign<0 else ap[:,i])
            v=(zi[:,None]+z[:,:i] if sign<0 else z[:,:i]-zi[:,None])
            for _ in range(7):
                if sign<0:
                    expo=beta-i
                    logf=expo*np.log1p(-2*r*zi)+np.sum(np.log1p(-r[:,None]*v),axis=1)
                    derivative=-2*expo*zi/(1-2*r*zi)-np.sum(v/(1-r[:,None]*v),axis=1)
                else:
                    expo=beta-n
                    logf=expo*np.log1p(2*r*zi)+np.sum(np.log1p(-r[:,None]*v),axis=1)
                    derivative=2*expo*zi/(1+2*r*zi)-np.sum(v/(1-r[:,None]*v),axis=1)
                residual=logf-np.log1p(-t)
                r-=residual/derivative
            assert np.all(r>0) and np.all(r<0.25)
            worst=max(worst,float(np.max(np.abs(residual))))
            answer.append(r)
    assert worst<2e-14
    return np.max(answer,axis=0),worst


def coefficients(n,beta):
    # Fractions keep every critical integer decision exact.
    c=min(beta,n-beta);sigma=2*c
    k=math.floor(2*beta) if beta<Fraction(n,2) else math.ceil(2*(n-beta))-1
    m=n-k;leading=Fraction(2**n,math.factorial(k))/sigma**(n-k)
    ratio=(Fraction(n,2)-(k+Fraction(m*(n+k+1),1)/sigma)/(n+1)
           if beta<Fraction(n,2) else Fraction(n,2)+(Fraction(m*(m+1),1)/sigma-k)/(n+1))
    return k,float(leading),float(leading*ratio)


def run_case(job):
    n,bstr,q,kind=job;bf=Fraction(bstr);beta=float(bf)
    z,w=directions(n,q);k,leading,second=coefficients(n,bf)
    sigma=2*min(beta,n-beta);v=z[:,k];y=z[:,:k]-v[:,None]
    a=z[:,:k].sum(axis=1)+(sigma-k)*v
    b=-.5*(y*y).sum(axis=1)
    if beta<n/2:
        b-=2*v*y.sum(axis=1)+sigma*v*v
    else:
        b+=sigma*v*v
    factor=2**n*math.factorial(n)
    leading_quad=factor/n*np.dot(w,a**(-n))
    second_quad=factor*np.dot(w,a**(-n)*(.5+b/a**2))
    rows=[]
    ts=[.01,.003,.001,.0001] if kind=='main' else [.01,.001]
    for t in ts:
        rho,residual=contact_roots(z,beta,t)
        scaled_volume=factor/n*np.dot(w,(rho/t)**n)
        # Common nodes cancel leading angular quadrature bias before t division.
        estimate=factor/n*np.dot(w,((rho/t)**n-a**(-n))/t)
        normalized_remainder=(estimate-second_quad)/t
        assert np.isfinite(normalized_remainder)
        rows.append({'t':t,'volume_over_tn':float(scaled_volume),
                     'second_estimate_same_nodes':float(estimate),
                     'remainder_over_t_nplus2_same_nodes':float(normalized_remainder),
                     'max_contact_log_root_residual':residual})
    return {'kind':kind,'n':n,'beta':bstr,'quadrature_order':q,'directions':len(w),
            'k':k,'L':leading,'T':second,'L_quadrature':float(leading_quad),
            'T_quadrature':float(second_quad),'rows':rows}


def main():
    started=time.perf_counter();base=Path(__file__).resolve().parents[1]
    if '--compare-only' in sys.argv:
        compare(base,json.loads((base/'results/higher_dimensional_python.json').read_text())['rows'])
        return
    table=json.loads((base/'results/higher_dimensional_wolfram.json').read_text())['rows']
    orders={1:1,2:128,3:64,4:32,5:20}
    jobs=[(r['n'],r['beta'],orders[r['n']],'main') for r in table]
    cross={3:['1/2','1','3/2','2','5/2'],4:['1','2','3'],5:['1','5/2','4']}
    jobs += [(n,b,{3:24,4:12,5:8}[n],'matlab-cross') for n,bs in cross.items() for b in bs]
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
        rows=list(pool.map(run_case,jobs))
    mainrows=[r for r in rows if r['kind']=='main']
    for row,ref in zip(mainrows,table):
        assert abs(row['L']/ref['L_numeric']-1)<2e-14
        assert abs(row['T']/ref['T_numeric']-1)<2e-14
    result={'kind':'non-certified-deterministic-numerical-evidence','numpy':np.__version__,
            'scipy':scipy.__version__,'precision':'IEEE754 double','random_seed':None,
            'root_iterations':7,'workers':8,'logical_cpus':os.cpu_count(),
            'elapsed_seconds':time.perf_counter()-started,
            'child_cpu_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime,
            'largest_child_peak_rss_bytes':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            'quadrature_note':'Tensor Gauss-Legendre on a Duffy map; no certified quadrature remainder.',
            'max_relative_L_quadrature_error':max(abs(r['L_quadrature']/r['L']-1) for r in mainrows),
            'max_relative_T_quadrature_error':max(abs(r['T_quadrature']/r['T']-1) for r in mainrows),
            'rows':rows}
    (base/'results/higher_dimensional_python.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
    compare(base,rows)


def compare(base,rows):
    matlab=base/'results/higher_dimensional_matlab.json'
    if matlab.exists():
        errors=[]
        for n,beta,q,t,value,*_ in json.loads(matlab.read_text())['rows']:
            ref=next(r for r in rows if r['kind']=='matlab-cross' and r['n']==n
                     and float(Fraction(r['beta']))==beta and r['quadrature_order']==q)
            target=next(r['volume_over_tn'] for r in ref['rows'] if r['t']==t)
            errors.append(abs(value/target-1))
        comparison={'kind':'cross-tool-comparison','cases':len(errors),
                    'max_relative_scaled_volume_difference':max(errors)}
        (base/'results/higher_dimensional_crosscheck.json').write_text(json.dumps(comparison,indent=2)+'\n')
        print(comparison)

if __name__=='__main__':
    main()
