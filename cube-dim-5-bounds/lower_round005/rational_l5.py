#!/usr/bin/env python3
"""One fully rational l<=5 witness certificate, inward from the locked C input.

The algebraic field is used only to obtain exact isolating lower bounds before
dyadic rounding.  Every box and every final volume below is a Fraction.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib, json, time, sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "lower"))
from certify_rectangular import field

D = 1 << 40
LIMIT_SECONDS, LIMIT_NODES = 120, 200000
started = time.monotonic()
nodes = 0

def dyadic_strict_lower(x, bounds):
    lo, _ = bounds(x.a)
    q = Q(lo.numerator * D // lo.denominator, D)
    assert q <= lo
    # Strictness is required even if an interval endpoint happens to be dyadic.
    if q == lo:
        q -= Q(1, D)
    assert q < lo
    return q, str(lo)

def inward_radius(n):
    """Largest dyadic m/D below (n/32)^(1/5), by integer arithmetic."""
    target = n * D**5
    lo, hi = 0, 1
    while 32 * hi**5 <= target:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if 32 * mid**5 <= target:
            lo = mid
        else:
            hi = mid
    if 32 * lo**5 == target:
        lo -= 1
    assert lo > 0 and 32 * lo**5 < target
    return Q(lo, D), lo

def union_volume(retained):
    global nodes
    @lru_cache(None)
    def rec(boxes):
        global nodes
        nodes += 1
        if nodes > LIMIT_NODES:
            raise RuntimeError("node cap")
        if time.monotonic() - started > LIMIT_SECONDS:
            raise TimeoutError("wall clock cap")
        if not boxes: return Q(0)
        if len(boxes[0]) == 2:
            out = Q(0); left = right = None
            for lo, hi in sorted(boxes):
                if left is None: left, right = lo, hi
                elif right < lo: out += right-left; left, right = lo, hi
                elif right < hi: right = hi
            return out + right-left
        cuts = sorted({x for b in boxes for x in b[:2]})
        total = Q(0)
        for lo, hi in zip(cuts, cuts[1:]):
            sub = tuple(sorted({b[2:] for b in boxes if b[0] <= lo and hi <= b[1]}))
            total += (hi-lo)*rec(sub)
        return total
    return rec(tuple(sorted(retained)))

def grid_volume(retained):
    segments=[]
    for j in range(5):
        cuts=sorted({b[2*j+k] for b in retained for k in (0,1)})
        ss=[]
        for lo,hi in zip(cuts,cuts[1:]):
            mask=0
            for i,b in enumerate(retained):
                if b[2*j] <= lo and hi <= b[2*j+1]: mask |= 1<<i
            if mask: ss.append((hi-lo,mask))
        segments.append(ss)
    @lru_cache(None)
    def rec(j,mask):
        if j==5:return Q(1)
        total=Q(0)
        for d,nxt in segments[j]:
            both=nxt if mask < 0 else mask & nxt
            if both: total += d*rec(j+1,both)
        return total
    return rec(0,-1), [len(x) for x in segments]

def main():
    F,u,v,bounds,_,_ = field()
    a1, a1_isolation = dyadic_strict_lower(1+u**3/2, bounds)
    tail, tail_isolation = dyadic_strict_lower((1+u)/2, bounds)
    spacing=(a1,tail,tail,tail,tail)
    raw=[]; radius_numerators={}
    for ls in product(range(1,6), repeat=5):
        n=prod(ls); r,m=inward_radius(n); radius_numerators[str(n)]=m
        box=[]
        for a,l in zip(spacing,ls):
            if l%2: lo,hi=Q(0),min(a/2,r-(l-1)*a/2)
            else: lo,hi=max(Q(0),l*a/2-r),a/2
            if hi <= lo: break
            box += [lo,hi]
        else: raw.append((ls,tuple(box),r,n))
    unique={b:(ls,r,n) for ls,b,r,n in raw}
    retained=[]
    for box,info in unique.items():
        if not any(other != box and all(other[2*j] <= box[2*j] and box[2*j+1] <= other[2*j+1] for j in range(5)) for other in unique):
            retained.append(box)
    # Directly re-check the strict radius inequality for every retained witness.
    for box in retained:
        ls,r,n=unique[box]
        assert 32*r**5 < n
    v1=union_volume(retained)
    v2, segment_counts=grid_volume(retained)
    assert v1 == v2
    ratio=32*v1
    c_lock=json.loads((ROOT/"lower_round003/lock_c_a2_equals_a3.json").read_text())
    c_upper=Q(c_lock["exact_interval"][1])
    out={
      "status":"PASS_exact_rational_l5" if ratio>c_upper else "NO_IMPROVEMENT_exact_rational_l5",
      "dimension":5,"length_family":"{1,...,5}^5","input_blocks":3125,
      "spacing":[str(x) for x in spacing],"spacing_denominator":str(D),
      "spacing_provenance":{"field":"u^5=2, v^5=3", "a1_lower_isolation":a1_isolation,
        "tail_lower_isolation":tail_isolation,"rounding":"strict floor to dyadic denominator 2^40"},
      "radius_rule":"maximal m/D with 32*m^5 <= n*D^5, decrement equality", "radius_numerators_by_product":radius_numerators,
      "raw_feasible_boxes":len(raw),"raw_feasible_boxes_with_length_5":sum(5 in ls for ls,_,_,_ in raw),
      "retained_boxes":len(retained),"retained_boxes_with_length_5":sum(5 in unique[b][0] for b in retained),
      "strict_radius_verified_for_all_retained":True,
      "union_engines_agree":True,"recursive_nodes":nodes,"coordinate_segment_counts":segment_counts,
      "ratio_fraction":str(ratio),"ratio_display_only":float(ratio),"C_exact_upper":str(c_upper),
      "strictly_above_C_exact_upper":ratio>c_upper,"wall_seconds":time.monotonic()-started,
      "limits":{"wall_seconds":LIMIT_SECONDS,"recursive_nodes":LIMIT_NODES},
      "collar_status":"not used; ±3 deferred to transfer audit", 
      "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }
    path=HERE/"rational_l5_results.json"; path.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__": main()
