#!/usr/bin/env python3
"""Independent exact recheck of the n=6 periodic-box lower-bound witness.

This does not import the original certificate implementation. It works in
Q(u,v), u^6=2, v^6=3, and orders elements by rational isolating intervals.
The volume sweep is a fresh recursive union computation.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib, json, time, math
from flint import arb, ctx

D = 10**45
n = 6

def floor_root(N, k):
    lo, hi = 0, 1 << ((N.bit_length() + k - 1)//k)
    while lo + 1 < hi:
        m = (lo + hi)//2
        if m**k <= N: lo = m
        else: hi = m
    return lo

ur = floor_root(2*D**n, n); vr = floor_root(3*D**n, n)
ul, uh, vl, vh = Q(ur,D), Q(ur+1,D), Q(vr,D), Q(vr+1,D)

class E:
    __slots__ = ('c',)
    def __init__(self, x=0, c=None):
        self.c = tuple(c) if c is not None else (Q(x),) + (Q(0),)*35
    def __hash__(self): return hash(self.c)
    def __eq__(self, y): return isinstance(y,E) and self.c == y.c
    def __add__(self, y):
        y = lift(y); return E(c=[a+b for a,b in zip(self.c,y.c)])
    __radd__ = __add__
    def __neg__(self): return E(c=[-a for a in self.c])
    def __sub__(self,y): return self + -lift(y)
    def __rsub__(self,y): return lift(y) + -self
    def __mul__(self,y):
        y=lift(y); out=[Q(0)]*36
        for i,a in enumerate(self.c):
            if not a: continue
            ai,aj=divmod(i,6)
            for j,b in enumerate(y.c):
                if not b: continue
                bi,bj=divmod(j,6); ii, jj=ai+bi, aj+bj
                out[(ii%6)*6+jj%6] += a*b*(2**(ii//6))*(3**(jj//6))
        return E(c=out)
    __rmul__=__mul__
    def __truediv__(self,k): return E(c=[a/Q(k) for a in self.c])
    def __pow__(self,k):
        z=E(1); b=self
        while k:
            if k&1:z=z*b
            b=b*b;k//=2
        return z
    @lru_cache(None)
    def interval(self):
        lo=Q(0);hi=Q(0)
        for k,a in enumerate(self.c):
            if not a:continue
            i,j=divmod(k,6); xlo=ul**i*vl**j; xhi=uh**i*vh**j
            if a>0:lo+=a*xlo;hi+=a*xhi
            else:lo+=a*xhi;hi+=a*xlo
        return lo,hi
    def __lt__(self,y):
        difference=lift(y)-self
        if not any(difference.c): return False
        lo,hi=difference.interval()
        if lo>0:return True
        if hi<0:return False
        # Some distinct high-degree endpoints have severe cancellation.  This
        # certified Arb fallback only decides their sign; all volume algebra
        # remains exact rational coefficient arithmetic.
        for prec in (1024,4096,16384):
            ctx.prec=prec; au=arb(2)**(arb(1)/6); av=arb(3)**(arb(1)/6); z=arb(0)
            for k,c in enumerate(difference.c):
                if c:
                    i,j=divmod(k,6)
                    z += arb(c.numerator)/arb(c.denominator)*au**i*av**j
            if 0 not in z: return z>0
        raise ArithmeticError('sign unresolved by exact interval or Arb')
    def __le__(self,y): return self==lift(y) or self<y

def lift(x): return x if isinstance(x,E) else E(x)
u=E(c=[Q(0)]*6+[Q(1)]+[Q(0)]*29)
v=E(c=[Q(0),Q(1)]+[Q(0)]*34)
assert u**6==E(2) and v**6==E(3)

def radical(N):
    a=b=0
    while N%2==0:a+=1;N//=2
    while N%3==0:b+=1;N//=3
    assert N==1
    return u**a*v**b

def maxe(a,b): return a if a==b else (b if a<b else a)
def mine(a,b): return a if a==b else (a if a<b else b)

def union_volume(boxes):
    @lru_cache(None)
    def rec(B):
        if not B:return E(0)
        if len(B[0])==2:
            ordered=sorted(B); a,b=ordered[0]; out=E(0)
            for c,d in ordered[1:]:
                if b<c:out+=b-a;a,b=c,d
                else:b=maxe(b,d)
            return out+b-a
        cuts=sorted({x for b in B for x in b[:2]}); total=E(0)
        for x,y in zip(cuts,cuts[1:]):
            active=tuple(sorted({b[2:] for b in B if b[0]<=x and y<=b[1]}))
            total+=(y-x)*rec(active)
        return total
    return rec(tuple(sorted(boxes)))

start=time.time()
spacing=[E(1)+u**4/2, (E(2)+u**5*v)/4]+[(E(1)+u)/2]*4
records=[]
for ell in product(range(1,5),repeat=6):
    r=radical(math.prod(ell))/2; box=[]
    for aa,k in zip(spacing,ell):
        lo,hi=(E(0), mine(aa/2,r-(k-1)*aa/2)) if k%2 else (maxe(E(0),k*aa/2-r),aa/2)
        if hi<=lo:break
        box += [lo,hi]
    else:
        starts=[-(k//2) if k%2 else 1-k//2 for k in ell]
        for j,(aa,k,s) in enumerate(zip(spacing,ell,starts)):
            for x in box[2*j:2*j+2]:
                for z in (s*aa,(s+k-1)*aa): assert -r <= x-z <= r
        records.append((ell,tuple(box)))

unique={b:e for e,b in records}; kept=[]
for b,e in unique.items():
    if any(c!=b and all(c[2*j]<=b[2*j] and b[2*j+1]<=c[2*j+1] for j in range(6)) for c in unique):continue
    kept.append((e,b))
V=union_volume(tuple(b for _,b in kept)); R=64*V
lo,hi=R.interval()
project=Path(__file__).resolve().parents[2]
ledger=json.loads((project/'results/bounds_ledger.json').read_text())
r4=ledger['inherited_lower']['interval'][1]
r4q=Q(r4)
assert E(c=[r4q]+[Q(0)]*35) < R
five_ledger=json.loads((project/'rounds/001/n5_readonly_snapshot.json').read_text())
r5=five_ledger['interval'][1]
assert E(c=[Q(r5)]+[Q(0)]*35) < R
coefficients={f'{i},{j}':str(c) for i in range(6) for j,c in enumerate(R.c[6*i:6*i+6]) if c}
out={'status':'pass','dimension':6,'candidate_formula':['1+u^4/2','(2+u^5*v)/4']+['(1+u)/2']*4,'u':'positive sixth root of 2','v':'positive sixth root of 3','raw_nonempty_boxes':len(records),'maximal_boxes':len(kept),'witness_counts':[list(e) for e,_ in kept],'polynomial_coefficients':coefficients,'R6_interval':[str(lo),str(hi)],'strictly_exceeds_R4_exact_upper':True,'R4_exact_upper':r4,'strictly_exceeds_active_R5_exact_upper':True,'R5_exact_upper':r5,'source_sha256':{'algebraic_bounds.py':'8b3cb2dd042332e52ef9ed3360c9cc8fce210ea37f2112584eb9f5a8d0c31030','calibration_and_transfer.md':'7ec33f05308b887cfb9923590d6ea727c6cf0ac66e4897356da551d2f04dbfda'},'elapsed_seconds':time.time()-start,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('recheck_n6_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
