#!/usr/bin/env python3
"""Pre-registered 15-point discovery line: a2*(1+k/10000), k=2..16."""
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib,json,time

TOL=1e-13
def ratio(k):
    u=2**(1/6);v=3**(1/6)
    a=(1+u**4/2,((2+u**5*v)/4)*(1+k/10000),*((1+u)/2 for _ in range(4)))
    raw=[]
    for ell in product(range(1,5),repeat=6):
        r=prod(ell)**(1/6)/2;b=[]
        for aj,q in zip(a,ell):
            lo,hi=(0.,min(aj/2,r-(q-1)*aj/2)) if q%2 else (max(0.,q*aj/2-r),aj/2)
            if hi<=lo+TOL:break
            b.extend((lo,hi))
        else:raw.append((ell,tuple(b)))
    unique={b:ell for ell,b in raw};keep=[]
    for b,e in unique.items():
        if any(c!=b and all(c[2*j]<=b[2*j]+TOL and b[2*j+1]<=c[2*j+1]+TOL for j in range(6)) for c in unique):continue
        keep.append((e,b))
    @lru_cache(None)
    def union(B):
        if not B:return 0.
        if len(B[0])==2:
            S=sorted(B);x,y=S[0];z=0.
            for a,b in S[1:]:
                if a>y+TOL:z,x,y=z+y-x,a,b
                elif b>y:y=b
            return z+y-x
        z=0.;cuts=sorted({x for b in B for x in b[:2]})
        for x,y in zip(cuts,cuts[1:]):
            B2=tuple(sorted({b[2:] for b in B if b[0]<=x+TOL and y<=b[1]+TOL}))
            z+=(y-x)*union(B2)
        return z
    return {'k':k,'raw_boxes':len(raw),'kept_boxes':len(keep),'ratio_float':64*union(tuple(sorted(b for _,b in keep)))}
def main():
    t=time.monotonic();out=[ratio(k) for k in range(2,17)];ref=2.06516698104057
    for x in out:x['margin_over_k1_display_upper']=x['ratio_float']-ref
    w=max(out,key=lambda x:x['ratio_float'])
    result={'status':'completed','scope':'numeric discovery only; no C6 claim','family':'a2*(1+k/10000), k=2..16, with R6-001 witness rule',
            'frozen_k1_reference_is_full_exact_interval_not_display':True,'k1_display_upper_for_diagnostic':ref,
            'candidates':out,'winner':w,'promotion_required':w['margin_over_k1_display_upper']>1e-7,
            'elapsed_seconds':time.monotonic()-t,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('second_spacing_line_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
