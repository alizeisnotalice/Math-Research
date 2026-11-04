#!/usr/bin/env python3
"""Exact integration of the cubic functional for four rectangular corner holes.

Piecewise affine candidate envelopes are integrated using Fraction only.
No sampled search is promoted to a universal inequality.
"""
from fractions import Fraction as F
from pathlib import Path
import json, random, time, platform


def objective(k,c,x,u):
    s,q=x-u,1-x-u
    mass=[a*b for a,b in k]
    a,b=k[c^1]; d,e=k[c^2]; h,j=k[c^3]
    return u*u+u+x*(1-x)+sum(mass)/2-mass[c] - max(a-s,0)*b-d*max(e-q,0)-max(h-s,0)*max(j-q,0)


def exact_cubic(k):
    """k[c]=(inward horizontal,vertical extent); c bit0=x sign,bit1=y sign."""
    k=[tuple(map(F,v)) for v in k]; total=F(0)
    for c in range(4):
        for L,R,ma,mb in [(F(0),F(1,2),F(1),F(0)),(F(1,2),F(1),F(-1),F(1))]:
            candidates=[(-ma,-mb),(ma,mb),(F(1),-k[c^1][0]),(F(-1),1-k[c^2][1]),(F(1),-k[c^3][0]),(F(-1),1-k[c^3][1])]
            cuts={L,R}
            for a,b in candidates:
                # validity and all positive-part branch boundaries
                for A,B in [(a+ma,b+mb),(a-ma,b-mb),(a-1,k[c^1][0]+b),(a+1,k[c^2][1]-1+b),(a-1,k[c^3][0]+b),(a+1,k[c^3][1]-1+b)]:
                    if A and L<-B/A<R: cuts.add(-B/A)
            cuts=sorted(cuts)
            for l,r in zip(cuts,cuts[1:]):
                mid=(l+r)/2; lines=[]
                for a,b in candidates:
                    if abs(a*mid+b)>ma*mid+mb: continue
                    y0=objective(k,c,l,a*l+b); y1=objective(k,c,r,a*r+b)
                    A=(y1-y0)/(r-l);B=y0-A*l
                    assert objective(k,c,mid,a*mid+b)==A*mid+B
                    lines.append((A,B))
                splits={l,r}
                for A,B in lines:
                    for C,D in lines:
                        if A!=C and l<(D-B)/(A-C)<r: splits.add((D-B)/(A-C))
                splits=sorted(splits)
                for lo,hi in zip(splits,splits[1:]):
                    A,B=max(lines,key=lambda z:z[0]*(lo+hi)/2+z[1])
                    total+=A*(hi*hi-lo*lo)/2+B*(hi-lo)
    return total


def main():
    start=time.monotonic(); rng=random.Random(20261105)
    # Exact square-hole formula, including all branch junctions and degeneracies.
    square=[]
    for k in [F(0),F(1,4),F(1,2),F(3,5),F(2,3),F(4,5),F(1),F(3,2),F(2)]:
        expected=(2-4*k**3 if k<=F(1,2) else 2+4*(k-1)*k*k if k<=F(2,3) else 6-18*k+22*k*k-8*k**3 if k<=1 else -2+12*k-8*k*k)
        got=exact_cubic([(k,k)]*4); assert got==expected
        square.append({'k':str(k),'B':str(got)})
    samples=[]
    for i in range(128):
        k=[(F(rng.randrange(17),8),F(rng.randrange(17),8)) for _ in range(4)]
        b=exact_cubic(k)
        # This finite assertion is a regression, not an all-profile theorem.
        assert b<=2
        samples.append({'corners':[[str(a),str(b)] for a,b in k],'B':str(b)})
    # Nonnegative defect identity for arbitrary pixel densities: exact rational check.
    checked=0
    for _ in range(512):
        x=F(rng.randrange(33),32);m=min(x,1-x);u=m*F(rng.randrange(-8,9),8);s=x-u;q=1-x-u
        density=[F(rng.randrange(9),8) for _ in range(16)];M=F(0);Tx=F(0);Ty=F(0);Txy=F(0);lower=F(0)
        for i in range(4):
            for j in range(4):
                a,b=F(i,2),F(i+1,2);c,d=F(j,2),F(j+1,2);f=density[4*i+j]
                dx=max(b-max(a,s),0);dy=max(d-max(c,q),0)
                M+=f*(b-a)*(d-c);Tx+=f*dx*(d-c);Ty+=f*(b-a)*dy;Txy+=f*dx*dy
                lower+=f*max(min(b,s)-a,0)*max(min(d,q)-c,0)
        phi=u*u+u+x*(1-x)+M-Tx-Ty-Txy
        defect=2*(m*m-u*u)+s*q-lower+2*Txy
        assert 2*m-phi==defect and defect>=0;checked+=1
    out={'python':platform.python_version(),'arithmetic':'fractions.Fraction; no numerical tolerance','seed':20261105,'square_regressions':square,'independent_rectangle_regressions':samples,'pixel_density_identity_checks':checked,'seconds':time.monotonic()-start,'scope':'Exact finite regressions and integrals; general theorem is proved in the manuscript. Independent corners are exploratory, not a universal upper bound.'}
    target=Path(__file__).resolve().parents[1]/'results/general_corners_python.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print({'square_cases':len(square),'four_corner_cases':len(samples),'density_checks':checked,'seconds':out['seconds']})
if __name__=='__main__':main()
