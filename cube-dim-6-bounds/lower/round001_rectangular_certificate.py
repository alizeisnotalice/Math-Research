#!/usr/bin/env python3
"""Exact finite six-dimensional rectangular-lattice lower certificate.

The computation uses Q[u,v]/(u^6-2,v^6-3), represented in the monomial
basis.  Sign decisions use rational isolating intervals; it uses no floating
point arithmetic.  It is the six-dimensional specialization of the locked
reference construction, rewritten here as a portable verifier.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import hashlib, json
from pathlib import Path

D, N = 10**45, 6

def floor_root(x, k):
    lo, hi = 0, 1
    while hi**k <= x: hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**k <= x: lo = mid
        else: hi = mid
    return lo

RU, RV = floor_root(2*D**N, N), floor_root(3*D**N, N)
UL, VL, UH, VH = Q(RU,D), Q(RV,D), Q(RU+1,D), Q(RV+1,D)
SZ = N*N

@lru_cache(None)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
@lru_cache(None)
def neg(a): return tuple(-x for x in a)
@lru_cache(None)
def mul(a,b):
    z = [Q(0)]*SZ
    for i,x in enumerate(a):
        if not x: continue
        ii,ij = divmod(i,N)
        for j,y in enumerate(b):
            if not y: continue
            ji,jj = divmod(j,N)
            z[((ii+ji)%N)*N + ((ij+jj)%N)] += x*y*2**((ii+ji)//N)*3**((ij+jj)//N)
    return tuple(z)
@lru_cache(None)
def bounds(a):
    lo=hi=Q(0)
    for k,c in enumerate(a):
        if not c: continue
        i,j=divmod(k,N); x,y=UL**i*VL**j, UH**i*VH**j
        if c > 0: lo,hi=lo+c*x,hi+c*y
        else: lo,hi=lo+c*y,hi+c*x
    return lo,hi
def positive(a):
    if not any(a): return False
    lo,hi=bounds(a)
    if lo>0: return True
    if hi<0: return False
    raise ArithmeticError('isolating intervals insufficient')

class F:
    __slots__=('a',)
    def __init__(self,n=0,a=None): self.a=tuple(a) if a is not None else (Q(n),)+(Q(0),)*(SZ-1)
    def __hash__(self): return hash(self.a)
    def __eq__(self,b): return self.a==(b.a if isinstance(b,F) else F(b).a)
    def __add__(self,b): return F(a=add(self.a,(b if isinstance(b,F) else F(b)).a))
    __radd__=__add__
    def __neg__(self): return F(a=neg(self.a))
    def __sub__(self,b): return self+-b
    def __rsub__(self,b): return F(b)+-self
    def __mul__(self,b): return F(a=mul(self.a,(b if isinstance(b,F) else F(b)).a))
    __rmul__=__mul__
    def __truediv__(self,k): return F(a=tuple(x/Q(k) for x in self.a))
    def __pow__(self,k):
        z=F(1)
        for _ in range(k): z=z*self
        return z
    def __lt__(self,b): return positive(((b if isinstance(b,F) else F(b))-self).a)
    def __le__(self,b): return self==b or self<b

u=F(a=(Q(0),)*N+(Q(1),)+(Q(0),)*(SZ-N-1))
v=F(a=(Q(0),Q(1))+(Q(0),)*(SZ-2))
assert u**N==F(2) and v**N==F(3)
def root(k):
    a=b=0
    while k%2==0: a,k=a+1,k//2
    while k%3==0: b,k=b+1,k//3
    assert k==1
    return u**a*v**b

@lru_cache(None)
def box_union(B):
    if not B: return F(0)
    if len(B[0])==2:
        ordered=sorted(B); a,b=ordered[0]; total=F(0)
        for c,d in ordered[1:]:
            if b<c: total,a,b=total+b-a,c,d
            else: b=max(b,d)
        return total+b-a
    total=F(0); cuts=sorted({x for box in B for x in box[:2]})
    for a,b in zip(cuts,cuts[1:]):
        child=tuple(sorted({box[2:] for box in B if box[0]<=a and b<=box[1]}))
        total += (b-a)*box_union(child)
    return total

def main():
    spacing=(1+u**4/2,(2+u**5*v)/4,*((1+u)/2 for _ in range(4)))
    records=[]
    for count in product(range(1,5),repeat=N):
        mass=1
        for x in count: mass*=x
        radius=root(mass)/2; box=[]
        for a,ell in zip(spacing,count):
            lo,hi=(F(0),min(a/2,radius-(ell-1)*a/2)) if ell%2 else (max(F(0),ell*a/2-radius),a/2)
            if hi<=lo: break
            box += [lo,hi]
        else:
            starts=[-(ell//2) if ell%2 else 1-ell//2 for ell in count]
            for j,(a,ell,start) in enumerate(zip(spacing,count,starts)):
                for x in box[2*j:2*j+2]:
                    for atom in (start*a,(start+ell-1)*a): assert -radius<=x-atom<=radius
            records.append((count,tuple(box)))
    unique={box:count for count,box in records}
    kept=[(count,box) for box,count in unique.items() if not any(other!=box and all(other[2*j]<=box[2*j] and box[2*j+1]<=other[2*j+1] for j in range(N)) for other in unique)]
    result=(2**N)*box_union(tuple(sorted(box for _,box in kept)))
    lo,hi=bounds(result.a); scale=10**14
    output={
      'status':'strict_finite_periodic_witness_union; transfer_to_L1_requires cited finite-block smoothing proof',
      'dimension':N,'parameter_formula':['1+u^4/2','(2+u^5*v)/4','(1+u)/2 x4'],
      'relations':['u^6=2','v^6=3'],'root_denominator':D,'root_lower_numerators':[RU,RV],
      'candidate_box_count':len(records),'noncontained_box_count':len(kept),
      'witness_counts':[list(c) for c,_ in kept], 'lower':str(lo),'upper':str(hi),
      'decimal_enclosure':['%d/10^14'%(lo.numerator*scale//lo.denominator),'%d/10^14'%((hi.numerator*scale+hi.denominator-1)//hi.denominator)],
      'polynomial_coefficients':{f'{i},{j}':str(result.a[i*N+j]) for i in range(N) for j in range(N) if result.a[i*N+j]},
      'self_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }
    print(json.dumps(output,indent=2))
if __name__=='__main__': main()
