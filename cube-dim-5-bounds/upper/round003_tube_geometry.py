#!/usr/bin/env python3
"""Exact finite geometry for the n=5 S2 tube diagnostic."""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations, product
from collections import Counter
from functools import reduce
from math import gcd
import hashlib, json
from pathlib import Path

def det(matrix):
    a=[list(row) for row in matrix]; previous=1; sign=1; n=len(a)
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        d=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=a[i][j]*d-a[i][k]*a[k][j]
                assert z%previous==0;a[i][j]=z//previous
        previous=d
        for i in range(k+1,n):a[i][k]=0
    return sign*a[-1][-1]

def interval_for_pair(s,t):
    low,high=Fraction(0),Fraction(1)
    for a,b in zip(s,t):
        if a==b:
            if a<0:return None
        elif a>b: # -1 to +1: -1+2lambda in [0,1]
            low=max(low,Fraction(1,2))
        else: # +1 to -1
            high=min(high,Fraction(1,2))
    return None if low>high else (low,high)

def main():
    tail=tuple(product((-1,1),repeat=4)); V=tuple((1,)+s for s in tail)
    histogram=Counter(); inverse_sum=Fraction(0)
    for triple in combinations(range(16),3):
        complement=tuple(j for j in range(16) if j not in triple)
        d=max(abs(det(tuple(V[i] for i in triple+pair))) for pair in combinations(complement,2))
        histogram[d]+=1;inverse_sum+=Fraction(1,d)
    assert histogram==Counter({16:160,32:240,48:160})
    assert inverse_sum==Fraction(125,6)
    nondegenerate=Counter(); point_traces=set(); segment_endpoints=set()
    for i,j in combinations(range(16),2):
        inter=interval_for_pair(tail[i],tail[j])
        if inter is None:continue
        lo,hi=inter
        point=lambda q:tuple(q*a+(1-q)*b for a,b in zip(tail[i],tail[j]))
        if lo==hi:point_traces.add(point(lo))
        else:
            changed=sum(a!=b for a,b in zip(tail[i],tail[j]))
            nondegenerate[changed]+=1;segment_endpoints.update((point(lo),point(hi)))
    assert nondegenerate==Counter({1:4,2:6,3:4,4:1})
    assert point_traces<=segment_endpoints
    payload={
      "triple_best_complement_det_histogram":{str(k):v for k,v in sorted(histogram.items())},
      "triple_inverse_det_sum":str(inverse_sum),
      "largest_label_inverse_det_sum":"125/2",
      "positive_face_S2_nondegenerate_segment_count_by_changed_coordinates":{str(k):v for k,v in sorted(nondegenerate.items())},
      "positive_face_S2_total_segment_length":"6+6*sqrt(2)+4*sqrt(3)",
      "positive_face_tube_bound":"for 0<delta<1/2: mu({dist_infinity(y,S2)<delta}) <= 256*delta^3-272*delta^4",
      "positive_face_tube_regions":"for each nonempty disagreement set J: fixed coordinates lie in [1-2delta,1], and the J coordinates have range <=2delta; summing binom(4,k)(2delta)^(4-k)[k(2delta)^(k-1)-(k-1)(2delta)^k] gives 256delta^3-272delta^4",
      "far_S2_heat_constant":"614656000/pi^5",
      "far_S2_T_bound":"T_epsilon(1,y) <= (307328000/pi^5)*delta^(-4)*epsilon^2 for dist_infinity(y,S2)>=delta",
      "scope":"The tube bound uses Euclidean radius 2delta to cover an infinity-norm delta tube around the 15 exact trace segments. It is a finite geometry/analytic diagnostic, not a C5 bound."
    }
    raw=json.dumps(payload,sort_keys=True,indent=2).encode()+b"\n";Path(__file__).with_name("round003_tube_geometry_results.json").write_bytes(raw)
    print(raw.decode(),end="");print("sha256",hashlib.sha256(raw).hexdigest())

if __name__=="__main__":main()
