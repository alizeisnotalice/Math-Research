#!/usr/bin/env python3
"""Exact real-window symmetry/mass checks, then uniform numerical quadrature."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json,time
import hashlib,importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent
THETA=[F(0),F(1,8),F(1,4),F(3,8)]
D=[F(1,8),F(3,16),F(1,4),F(5,16)]
DELTA=F(1,20);EPS=F(2,25)

def overlap_periodic(lo,hi,left,right):
    return sum((max(F(0),min(hi,right+j)-max(lo,left+j))for j in range((lo-right).__floor__()-1,(hi-left).__ceil__()+2)),F(0))
def segments(b):
    return [(F(0),F(1),EPS)]+[(THETA[b]+sg*D[b]-DELTA/2,THETA[b]+sg*D[b]+DELTA/2,(1-EPS)/(2*DELTA))for sg in (-1,1)]
def class_window(b,x,L,label_theta=None):
    """Exact selected-class integrals, unnormalized by L."""
    theta=THETA[b] if label_theta is None else label_theta
    lo,hi=x-L/2,x+L/2;out=[F(0),F(0)]
    for left,right,h in segments(b):
        for j in range((lo-right).__floor__()-1,(hi-left).__ceil__()+2):
            aa=max(lo,left+j);cc=min(hi,right+j)
            if aa>=cc:continue
            plus=overlap_periodic(aa,cc,theta,theta+F(1,2))
            out[0]+=h*plus;out[1]+=h*(cc-aa-plus)
    return out
def poly(rows):
    d=[F(1)]
    for plus,minus in rows:
        nd=[F(0)]*(len(d)+1)
        for k,v in enumerate(d):nd[k]+=v*minus;nd[k+1]+=v*plus
        d=nd
    return d
def value(x,L,w,same=False):
    ans=F(0)
    for b in range(4):
        d=poly([class_window(b,xx,L,F(0)if same else None)for xx in x])
        ans+=sum(v*w[b][k]for k,v in enumerate(d))/L**len(x)
    return ans
def point_value(x,w):
    ans=F(0)
    for b in range(4):
        density=F(1)
        for xx in x:
            density*=sum(h for left,right,h in segments(b)for shift in range(-2,3)if left+shift<=xx%1<right+shift)
        k=sum(int((xx-THETA[b])%1<F(1,2))for xx in x)
        ans+=density*w[b][k]
    return ans

def main():
    start=time.time();kernelchecks=0;same_class_differences=0
    points=[F(j,17)for j in range(17)]
    lengths=[F(1,5),F(3,4),F(5,4),F(3,2),F(9,4)]
    masses=[]
    for b in range(4):
        own=class_window(b,F(1,2),F(1))
        common=class_window(b,F(1,2),F(1),F(0))
        assert own==[F(1,2)]*2 and common==[F(1,2)]*2
        masses.append(dict(table=b,own=list(map(str,own)),common=list(map(str,common))))
        for x,L in product(points,lengths):
            z=class_window(b,x,L);mirror=class_window(b,2*THETA[b]-x,L)
            assert z==mirror[::-1]
            assert sum(z)==sum(mirror)
            kernelchecks+=1
            if z!=class_window(b,x,L,F(0)):same_class_differences+=1
    exacttilts=0
    for b,L,t in product(range(4),lengths,(-3,-1,0,1,3)):
        nodes=[F(2,7),2*THETA[b]-F(2,7),F(4,11),2*THETA[b]-F(4,11)]
        zz=[class_window(b,xx,L)for xx in nodes]
        aa=[sum(row)/L for row in zz];pp=[row[0]/sum(row)for row in zz]
        weights=[v**(1+t)for v in aa];den=sum(weights)
        assert sum(w*p for w,p in zip(weights,pp))/den==F(1,2)
        tau=4*sum(w*(p-F(1,2))**2 for w,p in zip(weights,pp))/den
        if L>=1:assert tau<=F(1,(2*int(L)+1)**2)
        exacttilts+=1
    # Exact finite input mass for nonnegative column-simplex weights.
    rng=np.random.default_rng(440044);masschecks=[]
    for n in (1,2,3,5,8):
        w=[[F(0)]*(n+1)for _ in range(4)]
        for k in range(n+1):
            nums=[int(a)for a in rng.integers(1,10,size=4)]
            for b in range(4):w[b][k]=F(nums[b],sum(nums))
        mass=sum(F(comb(n,k),2**n)*w[b][k]for b in range(4)for k in range(n+1))
        assert mass==1
        masschecks.append(dict(n=n,mass=str(mass)))
    # Explicit multiple positive table weights at a point, despite column simplex.
    n=2;w=[[F(int(b==(k%4)))for k in range(n+1)]for b in range(4)]
    point=[F(11,20),F(4,5)]
    point_sum=sum(w[b][sum(int((xx-THETA[b])%1<F(1,2))for xx in point)]for b in range(4))
    assert point_sum==2
    examples=[]
    whole_reflection_checks=0
    for n in (1,2,3):
        ww=[[F((b+1)*(k+1)+1)for k in range(n+1)]for b in range(4)]
        for k in range(n+1):
            den=sum(ww[b][k]for b in range(4))
            for b in range(4):ww[b][k]/=den
        xx=[F(2*i+1,11)for i in range(n)]
        v=value(xx,F(3,4),ww);same=value(xx,F(3,4),ww,True)
        examples.append(dict(n=n,own_labels_value=str(v),common_labels_value=str(same),different=v!=same))
        reflected=[ww[0][::-1],ww[3][:],ww[2][:],ww[1][:]]
        assert sum(F(comb(n,k),2**n)*reflected[b][k]for b in range(4)for k in range(n+1))==1
        for j in range(7):
            point=[F(2*i+2*j+1,37)for i in range(n)]
            assert point_value([-v for v in point],ww)==point_value(point,reflected)
            whole_reflection_checks+=1
    # Actual uniform quadrature: paired grid makes reflection test independent
    # from arbitrary query-law invariance assumptions. Use periodized box kernel.
    M=16384;x=(np.arange(M)+.5)/M;quads=[];max_prior_error=0.;max_long_trace=0.;max_long_excess=-1.
    def numeric_overlap(xx,L,left,right):
        z=np.zeros_like(xx)
        for j in range(-3,4):z+=np.maximum(0,np.minimum(xx+L/2,right+j)-np.maximum(xx-L/2,left+j))
        return z
    for b,L,t in product(range(4),[.2,.75,1.25,1.5,2.25],[-8,-1,0,1,8]):
        # Exact class intersections converted to base one-period segments.
        ss=[]
        for left,right,h in segments(b):
            for shift in range(-2,3):
                aa=max(F(0),left+shift);cc=min(F(1),right+shift)
                if aa>=cc:continue
                cuts=sorted({aa,cc}|{v for j in range(-1,2)for v in (THETA[b]+j,THETA[b]+F(1,2)+j)if aa<v<cc})
                for u,v in zip(cuts,cuts[1:]):
                    plus=((u+v)/2-THETA[b])%1<F(1,2)
                    ss.append((float(u),float(v),float(h),plus))
        z=np.zeros((M,2))
        for left,right,h,plus in ss:z[:,0 if plus else 1]+=h*numeric_overlap(x,L,left,right)/L
        a=z.sum(axis=1);p=z[:,0]/a
        ww=np.exp((1+t)*np.log(a)-np.max((1+t)*np.log(a)))
        mu=float(np.dot(ww,p)/ww.sum());max_prior_error=max(max_prior_error,abs(mu-.5))
        tau=float(np.dot(ww,(p-mu)**2)/ww.sum()/(mu*(1-mu)))
        if L>=1:
            q=int(L);cap=1/(2*q+1)**2
            max_long_trace=max(max_long_trace,tau);max_long_excess=max(max_long_excess,tau-cap)
        quads.append(dict(table=b,L=L,t=t,mu_plus=mu))
    assert max_prior_error<2e-12
    # Read-only main-kernel comparison, including wrapped b2 box and integer L.
    source=HERE.parent/'round044_symmetry/probe.py'
    spec=importlib.util.spec_from_file_location('main44',source)
    root=importlib.util.module_from_spec(spec);spec.loader.exec_module(root)
    realx=np.array([[float(v)for v in row]for row in [(F(0),F(1,7),F(1,2)),(F(1,20),F(19,20),F(3,4))]])
    lengths=[.25,.6,1.,1.25,1.5,2.,2.9]
    mainz=root.kernels(realx,lengths)
    main_error=0.;main_checks=0
    for b,i,j,li in product(range(4),range(3),range(2),range(len(lengths))):
        exact=class_window(b,F.from_float(float(realx[j,i])),F(str(lengths[li])))
        for outk,k in ((0,1),(1,0)):
            main_error=max(main_error,abs(mainz[b,i,j,li,outk]-float(exact[k]/F(str(lengths[li])))))
            main_checks+=1
    assert main_error<1e-12
    out=dict(status='PASS',exact_kernel_reflection_checks=kernelchecks,
        exact_symmetric_query_tilt_checks=exacttilts,
        different_own_vs_common_kernel_cases=same_class_differences,source_class_masses=masses,
        exact_input_mass_checks=masschecks,point_weight_sum_example=str(point_sum),
        own_vs_common_input_examples=examples,uniform_quadrature_cases=len(quads),
        exact_whole_reflection_checks=whole_reflection_checks,
        quadrature_nodes=M,max_tilted_prior_error_from_half=max_prior_error,
        max_long_window_nontrivial_trace=max_long_trace,max_excess_over_K2_cap=max_long_excess,
        main_kernel_fraction_checks=main_checks,main_kernel_max_error=main_error,
        main_script_SHA256=hashlib.sha256(source.read_bytes()).hexdigest(),
        scope='Exact mass and mirror kernel identities; finite uniform quadrature tilt checks. No new C_n estimate.',
        seconds=time.time()-start)
    (HERE/'CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
