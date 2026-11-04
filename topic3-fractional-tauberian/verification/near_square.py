#!/usr/bin/env python3
"""Exact certificates for the uniform near-square cubic upper bound.

The quadratic certificate covers all 24 order simplices, both x endpoints,
and all 15 degree-two Bernstein coefficients. Coverage identities use
independent rational cell integration; finite regressions are not the proof.
"""
from fractions import Fraction as F
from itertools import permutations, combinations
from collections import Counter
from pathlib import Path
import random, json, platform


def flat_R(x,u):
    y=1-x
    wx=min(u[1],u[3])+min(u[0],u[2])-2*x
    wy=min(u[2],u[3])+min(u[0],u[1])-2*y
    return sum((-2*x+u[c^2]+min(u[c^1],u[c^3]))*(-2*y+u[c^1]+min(u[c^2],u[c^3])) for c in range(4))-2*wx*wy


def flat_loss_over_x(x,v):
    u=[x*z for z in v]
    return (8*x-2*sum(u)-flat_R(x,u))/x


def certificate():
    rows=[]; histogram=Counter()
    for order in permutations(range(4)):
        vertices=[]
        for k in range(5):
            v=[F(0)]*4
            for j in range(4):v[order[j]]=F(-1 if j<k else 1)
            vertices.append(v)
        for endpoint in [F(0),F(1,2)]:
            def q(v):
                return flat_loss_over_x(F(1,2),v) if endpoint else 2*flat_loss_over_x(F(1,4),v)-flat_loss_over_x(F(1,2),v)
            coeff=[]
            for i in range(5):
                for j in range(i,5):
                    b=q(vertices[i]) if i==j else 2*q([(a+b)/2 for a,b in zip(vertices[i],vertices[j])])-(q(vertices[i])+q(vertices[j]))/2
                    assert b>=0;coeff.append(str(b));histogram[str(b)]+=1
            rows.append({'order':order,'x_endpoint':str(endpoint),'bernstein_upper_triangle':coeff})
    assert sum(histogram.values())==720
    return rows,dict(sorted(histogram.items()))


def intersection(rects):
    return max(min(r[1] for r in rects)-max(r[0] for r in rects),0)*max(min(r[3] for r in rects)-max(r[2] for r in rects),0)


def regression():
    rng=random.Random(20261105);checks=256
    for _ in range(checks):
        t=F(1,rng.randrange(100,1001));x=F(rng.randrange(65),64);y=1-x
        r=[1+t*F(rng.randrange(-8,9),8) for _ in range(4)]
        u=[F(rng.randrange(-8,9),16) for _ in range(4)]
        rects=[]
        for c in range(4):
            X=(-1 if c&1 else 1)*t*r[c]*x;Y=(-1 if c&2 else 1)*t*r[c]*y;h=F(1,2)+t*u[c]
            rects.append((X-h,X+h,Y-h,Y+h))
        triples=[intersection([rects[c] for c in S]) for S in combinations(range(4),3)]
        full=intersection(rects);positive=F(0)
        xs=sorted(set(a for r0 in rects for a in r0[:2]));ys=sorted(set(a for r0 in rects for a in r0[2:]))
        for a,b in zip(xs,xs[1:]):
            for c,d in zip(ys,ys[1:]):
                X=(a+b)/2;Y=(c+d)/2;k=sum(l<X<h and bot<Y<top for l,h,bot,top in rects)
                positive+=max(k-2,0)*(b-a)*(d-c)
        assert positive==sum(triples)-2*full
        W=[]
        for S in list(combinations(range(4),3))+[tuple(range(4))]:
            wx=min(((-1 if c&1 else 1)*r[c]*x+u[c]) for c in S)-max(((-1 if c&1 else 1)*r[c]*x-u[c]) for c in S)
            wy=min(((-1 if c&2 else 1)*r[c]*y+u[c]) for c in S)-max(((-1 if c&2 else 1)*r[c]*y-u[c]) for c in S)
            W.append((wx,wy))
        H=sum(a+b for a,b in W[:4])-2*sum(W[4]);R=sum(a*b for a,b in W[:4])-2*W[4][0]*W[4][1]
        assert H<=2*sum(u)-sum(r)
        assert 2+positive==4+t*H+t*t*R
    return checks


def main():
    rows,histogram=certificate();checks=regression()
    out={'python':platform.python_version(),'arithmetic':'Fraction exact; no tolerance','seed':20261105,'bernstein_coefficients':720,'coefficient_histogram':histogram,'certificate_rows':rows,'rational_coverage_regressions':checks,'scope':'The 720-coefficient certificate proves the stated quadratic inequality via the manuscript simplex decomposition. The 256 coverage regressions independently check formulas; localization and continuum transfer are analytic.'}
    target=Path(__file__).resolve().parents[1]/'results/near_square_python.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print({'coefficients':720,'histogram':histogram,'coverage_regressions':checks})
if __name__=='__main__':main()
