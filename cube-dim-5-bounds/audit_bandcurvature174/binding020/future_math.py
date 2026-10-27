"""UNEXECUTED future174 math primitives. No H/log/gradient evaluator.

Weights use positive P5/P3,not173's difference powers;corners clear one
positive integer denominator before taking the minimum. Separate primary
authorization of frozen actual173 outputs is required before real calls.
"""
from fractions import Fraction as F
from math import comb,lcm
from itertools import product
def positive_weight(a,b):
    a=F(a);b=None if b is None else F(b);assert a>0 and (b is None or b>a)
    B=F(0) if b is None else 1/b;d=1/a-B
    P=lambda n:sum((F(comb(n,k),k+1)*B**(n-k)*d**k for k in range(n+1)),F(0))
    return P(5)/P(3)
def source_q(original64,donor73,times,oldq):
    assert len(times)==len(donor73)==73 and times[24]=='1/64' and times[32]=='1/32'
    gL,gU=map(F,original64);ends=[tuple(map(F,x)) for x in donor73];dL,dU=ends[24];assert max(gL,dL)<=min(gU,dU)
    lower=max(gL,dL);assert lower>0;K=F(0);terms=[]
    for j in range(24,73):
        mass=ends[j][1]-(ends[j+1][0] if j<72 else F(0));assert mass>0
        weight=positive_weight(times[j],times[j+1] if j<72 else None);K+=mass*weight;terms.append((mass,weight))
    band=K/(16*lower);split=F(256)*(F(7,10)-F(8,15)*ends[32][0]/gU);q=min(F(oldq),band,F(512,3),split)
    assert 0<q<=F(oldq)
    return {'K':K,'lower64':lower,'band':band,'split':split,'q':q,'terms':terms}
def closed_integer_margin(point,delta115,gamma,q,box_lower,box_upper):
    point=list(map(F,point));delta=F(delta115);gamma=[tuple(map(F,x)) for x in gamma];q=F(q);lo=list(map(F,box_lower));hi=list(map(F,box_upper))
    assert point[0]==1 and q>0 and len(gamma)==4 and all(a<=b for a,b in gamma)
    D=lcm(*(x.denominator for x in [*point,delta,q,*lo,*hi,*[x for r in gamma for x in r]]));integer=lambda x:int(x*D)
    P=list(map(integer,point[1:]));A=integer(delta);Q=integer(q);G=[tuple(map(integer,r)) for r in gamma];L=list(map(integer,lo));U=list(map(integer,hi));out=[]
    for bits in product((0,1),repeat=4):
        X=[(U[j] if bits[j] else L[j])-P[j] for j in range(4)]
        num=A*D*D+D*sum(min(a*x,b*x) for (a,b),x in zip(G,X))-Q*sum(x*x for x in X);out.append((num,bits))
    value,bits=min(out);return F(value,D**3),bits
