#!/usr/bin/env python3
"""Exact finite periodic witness calculation for one alternating first-axis mass."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib, json, sys, time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
D=1<<40
LIMIT_SECONDS, LIMIT_NODES=120,200000

def strict_radius(mass):
    # r=m/D, with 32 r^5 < mass, using no floating point decisions.
    target=mass.numerator*D**5
    denom=32*mass.denominator
    lo,hi=0,1
    while denom*hi**5 <= target: hi*=2
    while lo+1<hi:
        mid=(lo+hi)//2
        if denom*mid**5 <= target: lo=mid
        else: hi=mid
    if denom*lo**5 == target: lo-=1
    assert lo>0 and denom*lo**5 < target
    return Q(lo,D),lo

def run(theta):
    start=time.monotonic(); nodes=0
    old=json.loads((ROOT/'lower_round005/rational_l5_results.json').read_text())
    a=tuple(Q(x) for x in old['spacing'])
    raw=[]; radius_numerators={}
    # The stated finite start range covers the receiver first-coordinate cell.
    for ls in product(range(1,5),repeat=5):
        l1=ls[0]
        for j in range(-5,4):
            # Python's (-1)**negative_integer is a float; retain exact parity.
            sign=Q(-1 if j % 2 else 1)
            first_mass=Q(l1) if l1%2==0 else Q(l1)+theta*sign
            mass=first_mass*prod(ls[1:])
            r,m=strict_radius(mass); radius_numerators.setdefault(str(mass),m)
            lo=max(Q(0),(j+l1-1)*a[0]-r)
            hi=min(2*a[0],j*a[0]+r)
            if hi<=lo: continue
            box=[lo,hi]
            for aa,l in zip(a[1:],ls[1:]):
                if l%2: qlo,qhi=Q(0),min(aa/2,r-(l-1)*aa/2)
                else: qlo,qhi=max(Q(0),l*aa/2-r),aa/2
                if qhi<=qlo: break
                box += [qlo,qhi]
            else: raw.append((j,ls,tuple(box),mass,r))
    # Identical boxes can be joined; then discard boxes contained in another.
    unique={b:(j,ls,mass,r) for j,ls,b,mass,r in raw}
    retained=[]
    for b in unique:
        if not any(o!=b and all(o[2*k]<=b[2*k] and b[2*k+1]<=o[2*k+1] for k in range(5)) for o in unique): retained.append(b)
    for b in retained:
        _,_,mass,r=unique[b]
        assert 32*r**5 < mass
    def check():
        nonlocal nodes
        nodes+=1
        if nodes>LIMIT_NODES: raise RuntimeError('node cap')
        if time.monotonic()-start>LIMIT_SECONDS: raise TimeoutError('wall clock cap')
    @lru_cache(None)
    def recursive(boxes):
        check()
        if not boxes:return Q(0)
        if len(boxes[0])==2:
            out=Q(0); left=right=None
            for lo,hi in sorted(boxes):
                if left is None:left,right=lo,hi
                elif right<lo:out+=right-left;left,right=lo,hi
                elif right<hi:right=hi
            return out+right-left
        total=Q(0); cuts=sorted({x for b in boxes for x in b[:2]})
        for lo,hi in zip(cuts,cuts[1:]):
            sub=tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]}))
            total+=(hi-lo)*recursive(sub)
        return total
    def grid():
        segments=[]
        for k in range(5):
            cuts=sorted({b[2*k+d] for b in retained for d in (0,1)})
            one=[]
            for lo,hi in zip(cuts,cuts[1:]):
                mask=0
                for i,b in enumerate(retained):
                    if b[2*k]<=lo and hi<=b[2*k+1]:mask|=1<<i
                if mask:one.append((hi-lo,mask))
            segments.append(one)
        @lru_cache(None)
        def rec(k,mask):
            if k==5:return Q(1)
            total=Q(0)
            for d,nextmask in segments[k]:
                both=nextmask if mask<0 else mask&nextmask
                if both:total+=d*rec(k+1,both)
            return total
        return rec(0,-1),[len(s) for s in segments]
    volume=recursive(tuple(sorted(retained)))
    grid_volume,segments=grid()
    assert volume==grid_volume
    ratio=8*volume
    c=json.loads((ROOT/'lower_round003/lock_c_a2_equals_a3.json').read_text())
    c_upper=Q(c['exact_interval'][1])
    return {"status":"PASS_finite_periodic_candidate" if ratio>c_upper else "NO_IMPROVEMENT_finite_periodic_candidate",
      "theta":str(theta),"spacing":list(map(str,a)),"spacing_source":"lower_round005 exact isolating-bound dyadic input",
      "input_length_vectors":1024,"first_axis_starts":[-5,3],"raw_boxes":len(raw),"unique_boxes":len(unique),"retained_boxes":len(retained),
      "radius_numerators_by_mass":radius_numerators,"strict_radius_verified_all_retained":True,
      "receiver_domain":"[0,2*a1] x product_{k=2..5}[0,ak/2]","period_mass":"2", "ratio_normalization":"8*union_volume",
      "union_engines_agree":True,"coordinate_segment_counts":segments,"recursive_nodes":nodes,"wall_seconds":time.monotonic()-start,
      "ratio_fraction":str(ratio),"ratio_display_only":float(ratio),"C_exact_upper":str(c_upper),"strictly_above_C_exact_upper":ratio>c_upper,
      "limits":{"wall_seconds":LIMIT_SECONDS,"recursive_nodes":LIMIT_NODES},
      "transfer_scope":"finite periodic calculation only; collar/smoothing transfer not supplied", 
      "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

def main():
    if len(sys.argv)!=2: raise SystemExit('usage: alternating_first_axis.py theta')
    theta=Q(sys.argv[1]); assert theta in (Q(1,8),Q(1,4),Q(1,2))
    try: out=run(theta)
    except (TimeoutError,RuntimeError) as e:
        out={"status":"COMPUTE_LIMIT", "theta":str(theta),"reason":str(e),"limits":{"wall_seconds":LIMIT_SECONDS,"recursive_nodes":LIMIT_NODES}}
    path=HERE/f'alternating_theta_{theta.numerator}_{theta.denominator}_results.json'
    path.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__':main()
