#!/usr/bin/env python3
"""Endpoint-grid union coverage for five cap-five boxes, without Q(u,v,w)."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib, json, sys, time

HERE=Path(__file__).resolve(); sys.path.insert(0,str(HERE.parents[1]))
import round001_rectangular_certificate as A

D=10**60; NODE_CAP=10**6; WALL=20*60
TARGETS=((3,5,5,5,5,5),(3,5,4,5,5,5),(3,5,5,4,5,5),
         (3,5,5,5,4,5),(3,5,5,5,5,4))

class Unresolved(Exception): pass
class Budget(Exception): pass

def floor_root(n):
    lo,hi=0,1
    while hi**6<=n: hi*=2
    while hi-lo>1:
        m=(lo+hi)//2
        if m**6<=n:lo=m
        else:hi=m
    return lo
def rootI(n):
    z=floor_root(n*D**6); return Q(z,D),Q(z+1,D)
def add(I,J):return I[0]+J[0],I[1]+J[1]
def sub(I,J):return I[0]-J[1],I[1]-J[0]
def mul(I,J):
    vals=(I[0]*J[0],I[0]*J[1],I[1]*J[0],I[1]*J[1]);return min(vals),max(vals)
def scale(I,c): return (I[0]*c,I[1]*c) if c>=0 else (I[1]*c,I[0]*c)
def power(I,k): return I[0]**k,I[1]**k
def div(I,J):
    assert J[0]>0
    vals=(I[0]/J[0],I[0]/J[1],I[1]/J[0],I[1]/J[1]);return min(vals),max(vals)

UV=(rootI(2),rootI(3)); W=rootI(5)
def field_bounds(f):
    """Rational outer enclosure of exact Q[u,v] element at D=1e60."""
    lo=hi=Q(0)
    for k,c in enumerate(f.a):
        if not c: continue
        i,j=divmod(k,6); I=mul(power(UV[0],i),power(UV[1],j))
        if c>=0: lo,hi=lo+c*I[0],hi+c*I[1]
        else: lo,hi=lo+c*I[1],hi+c*I[0]
    return lo,hi
def mrI(ell):
    n=prod(ell); ex=[]
    for p in (2,3,5):
        e=0
        while n%p==0:e,n=e+1,n//p
        ex.append(e)
    assert n==1
    out=(Q(1),Q(1))
    for I,e in zip((UV[0],UV[1],W),ex):out=mul(out,power(I,e))
    return out

u,v,F=A.u,A.v,A.F
SP=(F(1)+u**4/2,(u**3+u**4*v**3)/4,*((F(1)+u)/2 for _ in range(4)))
SPI=(field_bounds(x) for x in SP); SPI=tuple(SPI)
ZERO=F(0)

def old_box(ell):
    m=prod(ell); r=A.root(m)/2; b=[]
    for a,k in zip(SP,ell):
        lo,hi=(ZERO,min(a/2,r-(k-1)*a/2)) if k%2 else (max(ZERO,k*a/2-r),a/2)
        if hi<=lo:return None
        b += [lo,hi]
    return tuple(b)
def new_box(ell):
    """Endpoints are 0 exactly or rational enclosures; no w field object."""
    rr=scale(mrI(ell),Q(1,2)); b=[]
    for aI,k in zip(SPI,ell):
        half=scale(aI,Q(1,2))
        if k%2:
            lo=ZERO; hi=(min(half[0],sub(rr,scale(aI,Q(k-1,2)))[0]),
                          min(half[1],sub(rr,scale(aI,Q(k-1,2)))[1]))
        else:
            lo=(max(Q(0),sub(scale(aI,Q(k,2)),rr)[0]),
                max(Q(0),sub(scale(aI,Q(k,2)),rr)[1])); hi=(Q(1),Q(1))
            # physical x endpoint is a/2, not normalized 1
            hi=half
        if isinstance(lo,tuple): loE=lo
        else: loE=lo
        if isinstance(hi,tuple): hiE=hi
        else: hiE=hi
        loI=(Q(0),Q(0)) if loE is ZERO else loE
        if hiE[0]<=loI[1]: raise Unresolved('new box nonempty interval not separated')
        b += [(ZERO if loE is ZERO else loI),hiE]
    return tuple(b)

def exact_equal(x,y): return isinstance(x,F) and isinstance(y,F) and x==y
def interval(x): return field_bounds(x) if isinstance(x,F) else x
def cmp(x,y):
    """-1 x<y, 0 exact equality, +1 x>y; never tolerance-merges."""
    if exact_equal(x,y):return 0
    if isinstance(x,F) and isinstance(y,F):return -1 if x<y else 1
    xl,xh=interval(x); yl,yh=interval(y)
    if xh<yl:return -1
    if xl>yh:return 1
    raise Unresolved('endpoint order interval not separated')

def sortcuts(cuts):
    out=[]
    for x in cuts:
        found=False
        for y in out:
            if exact_equal(x,y):found=True;break
        if not found:out.append(x)
    # Insertion sort gives controlled unresolved reporting for mixed endpoints.
    ans=[]
    for x in out:
        pos=0
        while pos<len(ans) and cmp(ans[pos],x)<0:pos+=1
        if pos<len(ans) and cmp(ans[pos],x)==0:continue
        ans.insert(pos,x)
    return ans

def leq(x,y): return cmp(x,y)<=0

def cover_one(new, old, started):
    nodes=0
    @lru_cache(None)
    def rec(active, dim):
        nonlocal nodes
        nodes+=1
        if nodes>NODE_CAP:raise Budget('node cap')
        if time.monotonic()-started>WALL:raise Budget('wall time')
        if dim==6:return bool(active)
        lo,hi=new[2*dim],new[2*dim+1]
        cuts=[lo,hi]
        for i in active:
            a,b=old[i][2*dim],old[i][2*dim+1]
            # Only exact old cuts strictly inside the new dimension range.
            if cmp(lo,a)<0 and cmp(a,hi)<0:cuts.append(a)
            if cmp(lo,b)<0 and cmp(b,hi)<0:cuts.append(b)
        cuts=sortcuts(tuple(cuts))
        for left,right in zip(cuts,cuts[1:]):
            if cmp(left,right)>=0:raise Unresolved('nonpositive sorted slab')
            nxt=[]
            for i in active:
                a,b=old[i][2*dim],old[i][2*dim+1]
                if leq(a,left) and leq(right,b):nxt.append(i)
            if not nxt or not rec(tuple(nxt),dim+1):return False
        return True
    # Only boxes that can meet every new coordinate can enter the grid.
    active=[]
    for i,b in enumerate(old):
        if all(cmp(b[2*j+1],new[2*j])>=0 and cmp(b[2*j],new[2*j+1])<=0 for j in range(6)):
            active.append(i)
    ans=rec(tuple(active),0)
    return ans,nodes,len(active)

def main():
    t=time.monotonic()
    old=tuple(b for e in product(range(1,5),repeat=6) if (b:=old_box(e)) is not None)
    out={'scope':'symbolic endpoint-grid Boolean coverage; no volume arithmetic','D':str(D),'node_cap':NODE_CAP,
         'old_box_count':len(old),'targets':[]}
    try:
        for ell in TARGETS:
            n=new_box(ell); ok,nodes,active=cover_one(n,old,t)
            out['targets'].append({'ell':list(ell),'status':'covered' if ok else 'not_covered',
                                   'nodes':nodes,'initial_active_boxes':active})
        out['status']='pass' if all(x['status']=='covered' for x in out['targets']) else 'not_covered'
    except (Unresolved,Budget) as exc:
        out['status']='unresolved' if isinstance(exc,Unresolved) else 'budget_stop';out['reason']=str(exc)
    out['elapsed_seconds']=time.monotonic()-t
    out['code_sha256']=hashlib.sha256(HERE.read_bytes()).hexdigest()
    HERE.with_name('endpoint_grid_cover_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
