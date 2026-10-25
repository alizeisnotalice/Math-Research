#!/usr/bin/env python3
"""Fraction certificate: normalized old-inner union covers five new-outers."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib, json, time

HERE = Path(__file__).resolve()
TARGETS = ((3,5,5,5,5,5), (3,5,4,5,5,5), (3,5,5,4,5,5),
           (3,5,5,5,4,5), (3,5,5,5,5,4))

def floor_root(x):
    lo, hi = 0, 1
    while hi**6 <= x: hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**6 <= x: lo = mid
        else: hi = mid
    return lo

def roots(D):
    def one(n):
        z=floor_root(n*D**6); return Q(z,D),Q(z+1,D)
    return one(2),one(3),one(5)

def mul(I,J): return I[0]*J[0],I[1]*J[1]
def power(I,k): return I[0]**k,I[1]**k
def add(I,J): return I[0]+J[0],I[1]+J[1]
def scale(I,c): return I[0]*c,I[1]*c
def divide(I,J):
    assert I[0]>=0 and J[0]>0
    return I[0]/J[1],I[1]/J[0]

def mroot(ell, uvw):
    n=prod(ell); exps=[]
    for p in (2,3,5):
        e=0
        while n%p==0: e,n=e+1,n//p
        exps.append(e)
    assert n==1
    out=(Q(1),Q(1))
    for I,e in zip(uvw,exps): out=mul(out,power(I,e))
    return out

def spacings(uvw):
    u,v,_=uvw
    return (
      add((Q(1),Q(1)),scale(power(u,4),Q(1,2))),
      scale(add(power(u,3),mul(power(u,4),power(v,3))),Q(1,4)),
      *[scale(add((Q(1),Q(1)),u),Q(1,2)) for _ in range(4)])

def norm_box(ell, uvw, inner):
    """Inner old box or outer new box; all arithmetic is rational interval."""
    mr=mroot(ell,uvw); out=[]
    for e,a in zip(ell,spacings(uvw)):
        qlo,qhi=divide(mr,a) # q = 2r/a
        if e%2:
            lo=Q(0); hi=min(Q(1),(qlo if inner else qhi)-e+1)
        else:
            lo=max(Q(0),e-(qlo if inner else qhi)); hi=Q(1)
        if hi<=lo: return None
        out += [lo,hi]
    return tuple(out)

def union(B):
    @lru_cache(None)
    def rec(X):
        if not X: return Q(0)
        if len(X[0])==2:
            S=sorted(X); lo,hi=S[0]; ans=Q(0)
            for c,d in S[1:]:
                if c>hi: ans,lo,hi=ans+hi-lo,c,d
                elif d>hi: hi=d
            return ans+hi-lo
        ans=Q(0); cuts=sorted({z for b in X for z in b[:2]})
        for lo,hi in zip(cuts,cuts[1:]):
            if hi>lo:
                child=tuple(sorted({b[2:] for b in X if b[0]<=lo and hi<=b[1]}))
                ans+=(hi-lo)*rec(child)
        return ans
    return rec(tuple(sorted(B)))

def trial(D):
    uvw=roots(D)
    old=tuple(sorted(set(b for e in product(range(1,5),repeat=6)
                         if (b:=norm_box(e,uvw,True)) is not None)))
    new=tuple(norm_box(e,uvw,False) for e in TARGETS)
    if any(b is None for b in new):
        return {'D':str(D),'status':'unresolved_degenerate_outer_box'}
    Uold=union(old); Uboth=union(tuple(sorted(set(old+new))))
    return {'D':str(D),'status':'pass' if Uboth==Uold else 'unresolved_nonzero_residual',
            'old_inner_box_count':len(old),'new_outer_box_count':len(new),
            'residual':str(Uboth-Uold),'old_inner_volume':str(Uold),'joint_volume':str(Uboth)}

def main():
    t=time.monotonic(); out={'scope':'pure Fraction inward-old/outward-new cover','targets':[list(x) for x in TARGETS],'trials':[]}
    for D in (10**30,10**60):
        q=trial(D); out['trials'].append(q)
        if q['status']=='pass': break
    out['status']=out['trials'][-1]['status']; out['elapsed_seconds']=time.monotonic()-t
    out['code_sha256']=hashlib.sha256(HERE.read_bytes()).hexdigest()
    HERE.with_name('rational_union_cover_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
