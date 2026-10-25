#!/usr/bin/env python3
"""Independent finite audit: bisection roots, symbolic quotient ring, signed intersections."""
import hashlib, itertools, json, math, time
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import sympy as S
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'baseline'/'original_certificates.json'
X,Y=S.symbols('X Y')
D=10**35
@lru_cache(None)
def ri(k,n):
    # Deliberately independent, plain integer bisection, with certified brackets.
    if round(k**(1/n))**n==k:
        z=round(k**(1/n)); return F(z),F(z)
    target=k*D**n; lo=0; hi=max(2,k)*D
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**n<=target: lo=mid
        else: hi=mid
    assert lo**n<target<hi**n
    return F(lo,D),F(hi,D)
def add(*vs): return (sum(v[0] for v in vs),sum(v[1] for v in vs))
def mul(a,b):
    vs=[x*y for x in a for y in b]; return min(vs),max(vs)
def sc(a,c): return mul(a,(F(c),F(c)))
def pw(a,p):
    r=(F(1),F(1))
    for _ in range(p): r=mul(r,a)
    return r

def run(n,cert):
    u,v=ri(2,n),ri(3,n)
    ai=[add((F(1),F(1)),sc(pw(u,n-2),F(1,2))),sc(add((F(2),F(2)),mul(pw(u,n-1),v)),F(1,4))]+[sc(add((F(1),F(1)),u),F(1,2))]*(n-2)
    # Here the minimum spacing is identified with a proved disjoint interval.
    aminidx=min(range(n),key=lambda j:ai[j][0])
    assert all(ai[aminidx][1]<z[0] or z==ai[aminidx] for z in ai)
    t=(F(1),F(1))
    for z in ai:
        t=mul(t,(1/z[1]+1/(8*ai[aminidx][1]),1/z[0]+1/(8*ai[aminidx][0])))
    assert t[1]<1
    empty=[]; survivors=[]; lower_margin=None
    for ls in itertools.product(range(1,9),repeat=n):
        if max(ls)<=4: continue
        sd=ri(math.prod(ls),n)
        margins=[add(sd,sc(z,-(l-1))) for l,z in zip(ls,ai)]
        if any(z[1]<0 for z in margins):
            gap=-min(z[1] for z in margins)
            lower_margin=gap if lower_margin is None else min(lower_margin,gap)
            empty.append(ls)
        else:
            assert all(z[0]>0 for z in margins)
            assert all(l%2 for l in ls)
            assert all(z[1]<1 for z in margins)  # coordinate endpoint=(s-(l-1)a)/2 <1/2
            survivors.append(ls)
    # Exact old union: use a symbolic quotient ring and signed intersections,
    # algorithmically separate from the inherited coordinate-sweep verifier.
    def norm(z):
        z=S.Poly(S.expand(z),X,Y,domain=S.QQ)
        # Reduction as an explicit Groebner normal form, not inherited field code.
        return G.reduce(z.as_expr())[1]
    G=S.groebner([X**n-2,Y**n-3],X,Y,domain=S.QQ)
    a=[1+X**(n-2)/2,(2+X**(n-1)*Y)/4]+[(1+X)/2]*(n-2)
    @lru_cache(None)
    def enc(z):
        out=(F(0),F(0))
        for (i,j),c in S.Poly(z,X,Y).terms():
            out=add(out,sc(mul(pw(u,i),pw(v,j)),F(int(c.p),int(c.q))))
        return out
    def less(z,w):
        q=norm(w-z)
        if q==0:return False
        low,high=enc(q)
        assert low>0 or high<0,(n,q,low,high)
        return low>0
    def rt(k):
        p2=p3=0
        while k%2==0:k//=2;p2+=1
        while k%3==0:k//=3;p3+=1
        assert k==1
        return norm(X**p2*Y**p3)
    def makebox(ls):
        s=rt(math.prod(ls)); B=[]
        for l,z in zip(ls,a):
            if l%2:
                low=S.Integer(0); high=norm(s/2-(l-1)*z/2)
                if less(z/2,high):high=norm(z/2)
            else:
                low=norm(l*z/2-s/2)
                if less(low,0):low=S.Integer(0)
                high=norm(z/2)
            if not less(low,high):return None
            B.extend((low,high))
        return tuple(B)
    boxes=[makebox(ls) for ls in cert['witness_counts']]
    assert all(b is not None for b in boxes)
    all_nonempty=[]
    for ls in itertools.product(range(1,5),repeat=n):
        b=makebox(ls)
        if b is None:continue
        # Every old count box must lie in one retained certificate box.
        assert any(all(not less(b[2*j],c[2*j]) and not less(c[2*j+1],b[2*j+1]) for j in range(n)) for c in boxes),ls
        all_nonempty.append(ls)
    # Intern endpoints by exact expression; compute ordered ranks using independent interval comparisons.
    from functools import cmp_to_key
    cuts=[sorted({b[2*j+k] for b in boxes for k in (0,1)},key=cmp_to_key(lambda x,y: -1 if less(x,y) else (1 if less(y,x) else 0))) for j in range(n)]
    ids=[{x:i for i,x in enumerate(c)} for c in cuts]
    iboxes=[tuple(ids[j][b[2*j+k]] for j in range(n) for k in (0,1)) for b in boxes]
    signed={}
    for b in iboxes:
        delta={b:1}
        for c,coef in list(signed.items()):
            z=tuple(v for j in range(n) for v in (max(b[2*j],c[2*j]),min(b[2*j+1],c[2*j+1])))
            if all(z[2*j]<z[2*j+1] for j in range(n)):
                delta[z]=delta.get(z,0)-coef
        for z,coef in delta.items():
            signed[z]=signed.get(z,0)+coef
            if signed[z]==0:del signed[z]
    volume=S.Integer(0)
    for b,c in signed.items():
        term=S.Integer(c*2**n)
        for j in range(n):term*=cuts[j][b[2*j+1]]-cuts[j][b[2*j]]
        volume+=norm(term)
    volume=norm(volume)
    expected=sum(S.Rational(q)*X**int(key.split(',')[0])*Y**int(key.split(',')[1]) for key,q in cert['coefficients'].items())
    assert norm(volume-expected)==0
    vlo,vhi=enc(volume)
    oldlo,oldhi=map(F,cert['lower_bound_interval'])
    # Different root precision gives different enclosures, but exact coefficient equality above is the comparison.
    return dict(dimension=n,root_denominator=str(D),spacing_intervals=[[str(x),str(y)] for x,y in ai],tail_product_interval=[str(x) for x in t],tail_decimal=[float(x) for x in t],new_tuples=len(empty)+len(survivors),strictly_empty=len(empty),survivors=survivors,smallest_empty_margin=str(lower_margin),old_box_count=len(boxes),all_le4_nonempty=len(all_nonempty),all_le4_contained_in_certificate=True,signed_intersection_terms=len(signed),exact_coefficient_equality=True,volume_coefficients={f'{i},{j}':str(c) for (i,j),c in S.Poly(volume,X,Y).terms()},volume_interval=[str(vlo),str(vhi)],volume_decimal=[float(vlo),float(vhi)],inherited_interval_overlap=max(vlo,oldlo)<min(vhi,oldhi))
if __name__=='__main__':
    t=time.time(); certs=json.loads(BASE.read_text()); out=[run(c['dimension'],c) for c in certs]
    result={'input_sha256':hashlib.sha256(BASE.read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'engine':'SymPy '+S.__version__+' Groebner quotient ring, Fraction bisection; signed-intersection volume','dimensions':out,'elapsed_seconds':time.time()-t}
    (HERE/'002-independent-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
