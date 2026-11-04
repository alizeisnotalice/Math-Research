"""Exact lower tests for two asymmetric singularities. Standard library only.
The open interval (-1,1) is an analytic strict subset witness, not a sampled maximal set.
All square roots are bounded by verified integer inequalities at 180 bits.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import json
P=Path(__file__).resolve().parent
BITS=180

def root(x):
    x=Q(x);assert x>=0
    scale=1<<BITS
    k=isqrt(x.numerator*scale*scale//x.denominator)
    assert k*k*x.denominator<=x.numerator*scale*scale<(k+1)**2*x.denominator
    if k*k*x.denominator==x.numerator*scale*scale:return Q(k,scale),Q(k,scale)
    return Q(k,scale),Q(k+1,scale)

def primitive(x,positive,negative):
    lo,hi=root(abs(x))
    return (2*positive*lo,2*positive*hi) if x>=0 else (-2*negative*hi,-2*negative*lo)

def mass(a,b,z,positive,negative):
    al,ah=primitive(a-z,positive,negative);bl,bh=primitive(b-z,positive,negative)
    return bl-ah,bh-al

def mixture_mass(a,b,delta):
    x=mass(a,b,delta,1,2);y=mass(a,b,-delta,2,1)
    return x[0]+y[0],x[1]+y[1]

def test(delta,k,p):
    t=Q(1,2**k);h=t/(1+2*t);eps=delta*t**p
    assert (2+Q(8,5))*h<1 and h/(1-Q(8,5)*h)<t # strict common centered witnesses for the hole
    whole=mixture_mass(delta-eps,delta+eps,delta)
    hole=mixture_mass(delta-Q(8,5)*eps*h,delta+Q(2,5)*eps*h,delta)
    cost=(whole[0]-hole[1],whole[1]-hole[0]);assert cost[0]>0
    tl,th=root(t)
    norm=(hole[0]/cost[1]/th,hole[1]/cost[0]/tl)
    assert 0<norm[0]<=norm[1]
    # Explicit upper error for all |x|<=eps around +delta when eps<=delta.
    own=mass(delta-eps,delta+eps,delta,1,2)
    assert eps<=delta
    dl,dh=root(delta)
    other_upper=4*eps/dl
    assert whole[1]-own[0]<=other_upper
    return dict(delta=str(delta),k=k,epsilon_power=p,
                normalized_lower=str(norm[0]),normalized_upper=str(norm[1]),
                approximate_midpoint=float(sum(norm)/2))

def main():
    results=[test(Q(1,2**j),k,p) for j in (0,5,10,20) for k in range(4,33,2) for p in (1,2,4)]
    selected=[r for r in results if r['k']==32 and r['epsilon_power']==2]
    assert all(Q(r['normalized_lower'])>Q(21,20) for r in selected)
    a,b=root(10);a/=3;b/=3
    assert a>Q(21,20)>1
    summary=dict(status='passed',cases=len(results),bits=BITS,
                 halo='analytic strict subset: scaled (-1,1), endpoints discarded; not a full halo computation',
                 roots='integer isqrt checked inequalities',
                 exact_coefficient_gap='(sqrt(10)/3)^2-1=1/9',
                 coefficient_lower=str(a),coefficient_upper=str(b),
                 explicit_threshold_tests=selected,
                 last_k_normalized_range=[min(x['approximate_midpoint'] for x in selected),max(x['approximate_midpoint'] for x in selected)])
    assert Q(10,9)-1==Q(1,9)
    (P/'local_certificate_summary.json').write_text(json.dumps(summary,indent=2))
    (P/'local_certificate_full.json').write_text(json.dumps(results,indent=2))
    print('PASS',len(results),'strict subset witnesses and 180-bit bounds; 4 finite tests exceed 21/20')
    print(summary['last_k_normalized_range'])
if __name__=='__main__':main()
