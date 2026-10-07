#!/usr/bin/env python3
"""Exact local grid certificates; no numerical full-E integration."""
from fractions import Fraction as F
from math import ceil
import json


def radius(j):
    return F(4,2**j)


def grid_capture(x,r,s,h):
    lo=max(-6*h,(x-r)//s+1)
    hi=min(6*h,ceil((x+r)/s)-1)
    return lo,hi,max(0,hi-lo+1)


def case(q,batch):
    h=2**q; s=F(1,16*h); w=s/2
    near=(12*h+1)*w; far=1-near
    eps=s/1000; eta=s/100000; e=eps+eta
    assert F(3,8)<near<=F(13,32) and far>F(1,2)
    assert near/4<F(1,8)<near/2 and near<F(7,16)<F(15,32)<F(1,2)
    assert radius(q+7)==s/2 and radius(q+8)==s/4
    assert e<s/12
    band_q=e/(s/6)
    assert F(3,2)/(1+band_q)>1 and F(3,2)/(1-band_q)<2
    assert F(1,2)+s/(4*(F(5,6)*s))==F(4,5)<1
    antichain=[]
    for t in range(-h,h+1):
        x=(F(t)+F(1,6))*s
        tr={}
        for j in range(1,q+9):
            lo,hi,count=grid_capture(x,radius(j),s,h)
            tr[j]=count*w/(2*radius(j))
            assert abs(10-x)>radius(j)+e
            # Formula plus closest lattice boundary checks all near blobs.
            for endpoint in (x-radius(j),x+radius(j)):
                a=endpoint/s
                nearest=min(abs(a-(a//1)),abs(a-(a//1+1)))*s
                assert nearest>e
            for xx in (x-eps,x+eps):
                assert grid_capture(xx,radius(j),s,h)==(lo,hi,count)
        assert next(j for j,g in tr.items() if g>F(1,8))==2
        assert next(j for j,g in tr.items() if g>F(1,2))==q+8
        assert next(j for j,g in tr.items() if g>F(29,64))==4
        assert all(tr[j]==F(1,2) for j in range(4,q+8)) and tr[q+8]==1
        lo,hi,count=grid_capture(x,radius(4),s,h)
        assert (lo,hi,count)==(t-4*h+1,t+4*h,8*h)
        assert lo<=0<=hi
        antichain.append((lo,hi))
    assert len(set(antichain))==2*h+1
    assert all(b-a+1==8*h for a,b in antichain)
    assert F(1,8)*F(1,2)*2*eps==s/8000
    assert (s/8000)/w==F(1,4000)
    # Constants in the OPTIONAL full smooth allocation proof (identical even kernels).
    assert w/(s-2*eta)<1
    assert w/(4*(s/4-eta))<1
    assert 1/(10-F(3,8)-2*eta)<1
    assert far>4*eta
    assert far/(4*(far/2-eta))<1
    assert far/2+eta<2 and far/4>F(1,8)
    smooth_cap=F(3,4)+6*eta/s
    assert smooth_cap==F(37503,50000)<1
    return dict(batch=batch,q=q,h=h,actual_positive_cells=2*h+1,
                source_zero_mass=str(w),root_congestion_lower_bound=2*h+1,
                shared_beta_window=['7/16','15/32'],J=2,K_shared=4,K_half=q+8,
                atomic_full_E_explicit_capacity_upper='3/8',
                smooth_full_E_explicit_capacity_upper=str(smooth_cap),
                smooth_box_capacity_upper='1/4000',
                periodic_multi_blob_average_upper=str(w/(s-2*eta)),
                far_single_blob_superlevel_halo_radius=str(far/2),
                far_halo_J='1')


def main():
    batches=((0,1,2),(3,4,5),(6,7,8))
    cases=[case(q,batch) for batch,qs in enumerate(batches) for q in qs]
    return dict(status='passed',arithmetic='Fraction',cases=cases,
                scope='Exact local full-band, first-prefix, capture and smoothing constants. Full-E capacity bounds use separate analytic constructions, not sampled E integration.')


if __name__=='__main__':
    print(json.dumps(main(),indent=2))
