#!/usr/bin/env python3
"""Exact source-tree identities and complete integrals on a genuine 1D witness region."""
from fractions import Fraction as F
from pathlib import Path
import random,json,hashlib,sys

def tree_cost(leaves,V):
    b=sum(leaves,F(0));pair=raw=weighted=F(0)
    current=leaves[:]
    while len(current)>1:
        parent=[]
        for left,right in zip(current[::2],current[1::2]):
            v=left+right;parent.append(v)
            raw+=min(left,right)/V
            if b:
                pair+=left*right/(V*b)
                weighted+=v*min(left,right)/(V*b)
        current=parent
    terminal=sum((z*z for z in leaves),F(0))
    assert pair==(b*b-terminal)/(2*V*b) if b else pair==0
    assert pair<=weighted<=2*pair
    if b:assert weighted<=b/V
    return pair,raw,weighted

def capture(a,b,x,r):return max(F(0),min(b,x+r)-max(a,x-r))
def costs_at(x,L):
    ell=F(1,2**L);r=F(1,8)
    leaves=[capture(j*ell,(j+1)*ell,x,r) for j in range(2**L)]
    return tree_cost(leaves,2*r)

def full_pair_integral(L):
    # Piecewise quadratic. Breakpoints are every leaf endpoint +/- r.
    ell=F(1,2**L);r=F(1,8);a=F(1,4);b=F(3,4)
    points=sorted({a,b,*[z for j in range(2**L+1) for z in (j*ell-r,j*ell+r) if a<z<b]})
    total=F(0)
    for l,u in zip(points,points[1:]):
        # Exact Simpson integration, since the conditional root mass is constant.
        total+=(u-l)*(costs_at(l,L)[0]+4*costs_at((l+u)/2,L)[0]+costs_at(u,L)[0])/6
    return total

def raw_layer_integral(k):
    ell=F(1,2**k);r=F(1,8);a=F(1,4);b=F(3,4)
    points={a,b}
    for j in range(2**k):
        l=j*ell;mid=l+ell/2;u=l+ell
        for y in (l,mid,u):
            for z in (y-r,y+r):
                if a<z<b:points.add(z)
    # On initial cells, each child capture is affine; min may add one root.
    cuts=sorted(points)
    for l,u in zip(cuts,cuts[1:]):
        for j in range(2**k):
            left=j*ell;mid=left+ell/2;right=left+ell
            def difference(x):return capture(left,mid,x,r)-capture(mid,right,x,r)
            dl,du=difference(l),difference(u)
            if dl*du<0:points.add(l-dl*(u-l)/(du-dl))
    def value(x):
        return sum((min(capture(j*ell,j*ell+ell/2,x,r),capture(j*ell+ell/2,(j+1)*ell,x,r)) for j in range(2**k)),F(0))/(2*r)
    cuts=sorted(points)
    return sum(((u-l)*(value(l)+value(u))/2 for l,u in zip(cuts,cuts[1:])),F(0))

def main():
    batches=[]
    for batch,(depth,k,L) in enumerate(((3,2,2),(6,4,4),(9,6,6))):
        rng=random.Random(710000+batch);checks=0
        for trial in range(8):
            weights=[rng.randint(0,10000) for _ in range(2**depth)]
            total=sum(weights);q=[F(w,total) for w in weights]
            captured=[v*F(rng.randint(0,16),16) for v in q]
            V=F(rng.randint(1,9),16)
            pair,raw,weighted=tree_cost(captured,V)
            if sum(captured):
                b=sum(captured);m=b-max(captured);d=b-sum(z*z for z in captured)/b
                assert m<=d<=2*m
                # Fine-group impurity is 2V times the tree pair cost.
                assert d==2*V*pair
            checks+=1
        assert tree_cost([F(0)]*(2**depth),F(1))==(F(0),F(0),F(0))
        singleton=[F(0)]*(2**depth);singleton[-1]=F(1)
        assert tree_cost(singleton,F(1))==(F(0),F(0),F(0))
        checks+=2
        # Balanced capture: the local-node denominator still pays every depth.
        equal=[F(1,2**depth)]*(2**depth)
        pair_bal,raw_bal,weighted_bal=tree_cost(equal,F(1))
        assert raw_bal==F(depth,2)
        assert weighted_bal==1-F(1,2**depth)
        assert pair_bal==weighted_bal/2
        pair_exact=full_pair_integral(L);ell=F(1,2**L)
        pair_formula=(1-4*ell+F(16,3)*ell*ell)/4
        assert pair_exact==pair_formula
        raw_exact=raw_layer_integral(k)
        assert raw_exact==F(1,4)-F(1,2**(k+1))
        large_depth=(16,64,256)[batch]
        # Closed-form tail of the proven per-level formula; no huge tree enumerated.
        raw_tail=F(large_depth-4,4)-F(1,16)+F(1,2**large_depth)
        assert raw_tail>=F(large_depth-4,8)
        terminal=F(1,2**large_depth)
        pair_large=(1-4*terminal+F(16,3)*terminal**2)/4
        assert 0<pair_large<F(1,4)
        # Exact fractional-label benchmark for the endpoint-equivalence audit.
        a=F(5,4);N=2**depth
        assert a-max([a/N]*N)==a*(1-F(1,N))
        # Spatially disjoint periodic labels on the same continuous plateau.
        period=F(1,16*2**batch);groups=2**(batch+1);w=F(1,4)
        label_data=[[] for _ in range(groups)]
        for j in range(int(1/period)):
            for i in range(groups):
                label_data[i].append((j*period+i*period/groups,j*period+(i+1)*period/groups))
        stripe_checks=0
        for x in (F(1,4),F(1,3),F(2,5),F(1,2),F(3,4)):
            values=[sum((capture(l,u,x,w/2) for l,u in intervals),F(0)) for intervals in label_data]
            assert values==[w/groups]*groups
            tp,tr,tw=tree_cost(values,w)
            assert tp==(1-F(1,groups))/2
            assert tr==F(batch+1,2)
            stripe_checks+=1
        batches.append(dict(batch=batch,seed=710000+batch,abstract_tree_depth=depth,tree_checks=checks,balanced_tree_checks=1,stripe_checks=stripe_checks,stripe_groups=groups,
            integrated_pair_depth=L,pair_integral=str(pair_exact),integrated_raw_layer=k,raw_integral=str(raw_exact),
            formula_only_depth=large_depth,raw_tail=str(raw_tail),pair_formula_below_quarter=True))
    p=Path(__file__)
    print(json.dumps(dict(round=71,status='passed',python=sys.version.split()[0],batches=batches,
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()},
        scope='Abstract tree algebra plus exact 1D spatial integrals; large depths use proved closed formulas, not enumerated trees; no endpoint bound proved.',
        main_theorem_proved=False),indent=2))
if __name__=='__main__':main()
