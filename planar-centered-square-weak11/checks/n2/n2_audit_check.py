from fractions import Fraction as F
from pathlib import Path
import hashlib
import json


def add(p,q):
    z=p.copy()
    for ij,c in q.items():
        z[ij]=z.get(ij,F(0))+c
    return {ij:c for ij,c in z.items() if c}


def scale(p,c):
    return {ij:v*c for ij,v in p.items() if v*c}


def sub(p,q):
    return add(p,scale(q,F(-1)))


def mul(p,q):
    z={}
    for (i,j),c in p.items():
        for (k,l),d in q.items():
            z[(i+k,j+l)]=z.get((i+k,j+l),F(0))+c*d
    return {ij:c for ij,c in z.items() if c}


def rem_circle(p):
    p=p.copy()
    while any(i>=2 for i,j in p):
        ij=next(ij for ij in p if ij[0]>=2)
        c=p.pop(ij);i,j=ij
        p=add(p,{(i-2,j):c,(i-2,j+2):-c})
    return p


one={(0,0):F(1)}
zero={}
aa={(1,0):F(1)}
bb={(0,1):F(1)}
b2=mul(bb,bb)
ts=[zero,scale(sub(one,aa),F(1,2)),scale(sub(one,bb),F(1,2)),
    scale(add(one,bb),F(1,2)),scale(add(one,aa),F(1,2)),one]
la=[aa,aa,scale(add(aa,bb),F(1,2)),scale(sub(aa,bb),F(1,2)),zero,zero]
lb=[bb,bb,bb,zero,zero,zero]
HH=scale(mul(bb,sub(add(aa,bb),one)),F(1,2))
H0=scale(sub(one,aa),F(1,2))
expected=[
    [zero],
    [scale(H0,-1),sub(scale(H0,-1),scale(b2,F(1,4)))],
    [scale(HH,-1),
     scale(add(mul(sub(one,aa),sub(one,bb)),scale(b2,2)),F(-1,4)),
     sub(scale(HH,-1),scale(b2,F(1,4)))],
    [HH,sub(HH,scale(mul(sub(one,aa),sub(one,bb)),F(1,4))),
     scale(b2,F(1,4)),sub(HH,scale(b2,F(1,4)))],
    [H0,scale(b2,F(1,4)),scale(mul(H0,add(one,bb)),F(1,2)),
     scale(mul(H0,sub(one,bb)),F(1,2)),sub(H0,scale(b2,F(1,4)))],
    [zero]*6]
table=[]
for i in range(6):
    for j in range(i+1):
        actual=sub(sub(mul(sub(one,ts[i]),sub(one,ts[j])),
                       mul(la[i],la[j])),mul(lb[i],lb[j]))
        residual=rem_circle(sub(actual,expected[i][j]))
        assert not residual,(i,j,residual)
        table.append([i,j,True])
assert not rem_circle(sub(sub(HH,H0),
                         scale(mul(sub(one,aa),sub(aa,bb)),F(1,2))))


def union_area(rects):
    rects=[r for r in rects if r[0]<r[1] and r[2]<r[3]]
    if not rects:
        return F(0)
    xx=sorted({x for l,r,b,t in rects for x in (l,r)})
    yy=sorted({y for l,r,b,t in rects for y in (b,t)})
    out=F(0)
    for l,r in zip(xx,xx[1:]):
        for b,t in zip(yy,yy[1:]):
            x,y=(l+r)/2,(b+t)/2
            if any(ll<x<rr and bb<y<tt for ll,rr,bb,tt in rects):
                out+=(r-l)*(t-b)
    return out


def ell(c,z):
    return min(c,max(F(0),(1+c)/2-z))


geometry=0
for z in [F(1,100),F(1,10),F(1,5),F(2,5)]:
    a=(1-z*z)/(1+z*z);b=2*z/(1+z*z)
    assert a>=b>0 and a*a+b*b==1
    knots=[F(0),(1-a)/2,(1-b)/2,(1+b)/2,(1+a)/2,F(1),F(3,2)]
    probes=sorted(set(knots+[(knots[i]+knots[i+1])/2
                            for i in range(len(knots)-1)]))
    for x in probes:
        for y in probes:
            r=[(-a/2,a/2,-a/2,a/2),
               (x-b/2,x+b/2,y-b/2,y+b/2),
               (x-F(1,2),F(1,2),y-F(1,2),F(1,2))]
            actual=union_area(r)
            formula=1+max(1-x,0)*max(1-y,0)-ell(a,x)*ell(a,y)-ell(b,x)*ell(b,y)
            assert actual==formula,(z,x,y,actual,formula)
            H=b*(a+b-1)/2
            assert actual<=1+H
            geometry+=1

lip=0
for big in [F(0),F(1,4),F(1),F(2)]:
    for small in [F(0),big/2,big]:
        for dx in [F(0),F(1,5),F(1),F(3)]:
            for dy in [F(0),F(1,3),F(1),F(3)]:
                loss=max(big-dx,0)*max(big-dy,0)-max(small-dx,0)*max(small-dy,0)
                assert 0<=loss<=big*big-small*small
                lip+=1
sharp=[]
for z in [F(1,10),F(1,100),F(1,1000)]:
    a=(1-z*z)/(1+z*z);b=2*z/(1+z*z)
    H=b*(a+b-1)/2
    val=1+H/(b*b)
    assert val==F(3,2)-b/(2*(1+a))
    sharp.append({'circle_parameter':str(z),'deletion_loss_over_mass':str(val)})

out={'scope':'independent finite algebra and rational geometric checks; analytic proof audited separately',
     'all_21_polynomial_identities_mod_circle':table,
     'rational_three_rectangle_union_checks':geometry,
     'box_mass_Lipschitz_checks':lip,
     'sharp_deletion_limit_examples':sharp,
     'all_assertions_passed':True,
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
