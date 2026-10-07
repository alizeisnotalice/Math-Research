#!/usr/bin/env python3
"""Exact layer-coupling audit. Synthetic marginals, NOT actual ball examples."""
from fractions import Fraction as F
from math import comb, factorial
import json, random, sys

def pos(x): return max(F(0), x)

def integrate(fun, lo, hi):
    return (hi-lo)*(fun(lo)+4*fun((lo+hi)/2)+fun(hi))/6

def cost(S,h,C,segments):
    independent=monotone=F(0)
    assert h>0 and segments[0][0]==0 and segments[-1][1]==1
    previous=F(0);last=F(0)
    for lo,hi,b0,b1 in segments:
        assert lo==previous and hi>lo and b1>=b0>=last>=0
        previous=hi;last=b1
        slope=(b1-b0)/(hi-lo)
        def b(u): return b0+slope*(u-lo)
        def product(u):
            a=C+b(u)
            return (pos(S+h-a)**2-pos(S-a)**2)/2
        def mono(u): return h*pos(S+h*u-C-b(u))
        cuts={lo,hi}
        if slope:
            for target in (S-C,S+h-C):
                root=lo+(target-b0)/slope
                if lo<root<hi:cuts.add(root)
        if h!=slope:
            root=(C+b0-slope*lo-S)/(h-slope)
            if lo<root<hi:cuts.add(root)
        cuts=sorted(cuts)
        for a,z in zip(cuts,cuts[1:]):
            independent+=integrate(product,a,z)
            monotone+=integrate(mono,a,z)
    # Separate tail-excess computation; target law may have atoms and jumps.
    def target_tail(t):
        total=F(0)
        for lo,hi,b0,b1 in segments:
            if b0==b1:total+=(hi-lo)*int(C+b0>t)
            elif t<=C+b0:total+=hi-lo
            elif t<C+b1:total+=(hi-lo)*(C+b1-t)/(b1-b0)
        return h*total
    def difference(t): return min(h,pos(S+h-t))-target_tail(t)
    knots=sorted({F(0),S,S+h}|{C+b for segment in segments for b in segment[2:]})
    tail=F(0)
    for a,z in zip(knots,knots[1:]):
        # Interior points avoid one-sided ambiguity at target atoms.
        p=(2*a+z)/3;q=(a+2*z)/3
        slope=(difference(q)-difference(p))/(q-p)
        def affine(t):return difference(p)+slope*(t-p)
        cuts=[a,z]
        if slope:
            root=p-difference(p)/slope
            if a<root<z:cuts.insert(1,root)
        for left,right in zip(cuts,cuts[1:]):
            tail+=(right-left)*pos(affine((left+right)/2))
    assert tail==monotone
    gain=independent-monotone
    assert F(0)<=gain<=h*h/6
    return dict(product=str(independent),optimal=str(monotone),tail_formula=str(tail),gain=str(gain),gain_over_h_squared=str(gain/(h*h)))

def cases():
    for S,h,C in [(F(0),F(1),F(0)),(F(2),F(1,4),F(1)),(F(100),F(1,1024),F(3,2))]:
        yield 0,'sharp_uniform',S,h,C,[(F(0),F(1),S-C,S-C+h)]
    for target in (F(0),F(1,8),F(1,2),F(1),F(3)):
        yield 0,'constant',F(1,4),F(1,2),F(1,8),[(F(0),F(1),target,target)]
    for width in (2,4,8,16):
        seg=[(F(i,width),F(i+1,width),F(i,width),F(i,width)) for i in range(width)]
        yield 1,'staircase',F(0),F(1),F(0),seg
    for offset in (F(0),F(1,2),F(10)):
        yield 1,'atoms_and_ramps',offset,F(1,2),F(1,4),[(F(0),F(1,4),F(0),F(0)),(F(1,4),F(3,4),F(1,8),F(2)),(F(3,4),F(1),F(3),F(3))]
    rng=random.Random(450710)
    for count in (8,16,32):
        for repeat in range(5):
            vals=sorted(F(rng.randrange(129),32) for _ in range(2*count))
            segments=[(F(i,count),F(i+1,count),vals[2*i],vals[2*i+1]) for i in range(count)]
            yield 2,'random_monotone',F(rng.randrange(65),8),F(1,2**repeat),F(rng.randrange(9),8),segments

def ball_benchmarks():
    """Exact lens upper bounds for the actual separated two-source exit shells.

    These certify bounds, NOT exact integration of the high-dimensional deficit.
    Unit-ball center separation 1/2 is a relaxation of true |z-a_i|>2^(-1/n).
    """
    batches=((1,3,5),(9,17,33),(65,129,257,513,1025))
    rows=[]
    for batch,ns in enumerate(batches):
        for n in ns:
            m=(n-1)//2
            beta=F(factorial(2*m+1),2**(2*m+1)*factorial(m)**2)
            loss=2*beta*sum(((-1)**i*comb(m,i)*F(1,4)**(2*i+1)/F(2*i+1) for i in range(m+1)),F(0))
            eps=1-loss
            assert 0<eps<=min(F(1),F(16,n+2))
            limits=[]
            for C in (F(0),F(1,4),F(49,100),F(1,2)):
                lower=pos(F(1,2)-C-eps)**2
                upper=pos(F(1,2)-C)**2
                assert 0<=lower<=upper
                limits.append(dict(C=str(C),D_over_X_lower=str(lower),D_over_X_upper=str(upper),lower_display=float(lower)))
            rows.append(dict(batch=batch,n=n,epsilon=str(eps),epsilon_display=float(eps),markov_bound=str(min(F(1),F(16,n+2))),bounds=limits))
    assert F(rows[-1]['bounds'][2]['D_over_X_lower'])>0
    return dict(scope='exact cap upper certificate for analytically derived actual two-atom band; not numerical integration of D',case_count=len(rows),batch_counts=[len(ns) for ns in batches],cases=rows)

def main():
    rows=[]
    for batch,label,S,h,C,segments in cases():
        row=cost(S,h,C,segments)
        if label=='sharp_uniform':assert F(row['gain_over_h_squared'])==F(1,6) and F(row['optimal'])==0
        if label=='constant':assert F(row['gain'])==0
        rows.append(dict(batch=batch,label=label,S=str(S),h=str(h),C=str(C),segment_count=len(segments),**row))
    result=dict(status='passed',arithmetic='Fraction exact',scope='synthetic one-layer marginals, not Euclidean realizations',case_count=len(rows),batch_counts=[sum(c['batch']==b for c in rows) for b in range(3)],max_gain_over_h_squared=str(max(F(c['gain_over_h_squared']) for c in rows)),seed=450710,cases=rows,ball_benchmarks=ball_benchmarks())
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
