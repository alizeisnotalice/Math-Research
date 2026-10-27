"""Independent172 arithmetic. Import is inert; real calls need primary approval."""
from fractions import Fraction as F
from math import comb,lcm
from itertools import product
def weight(a,b):
    a=F(a);b=None if b is None else F(b);assert a>0 and (b is None or b>a)
    B=F(0) if b is None else 1/b;d=1/a-B
    P=lambda n:sum((F(comb(n,k),k+1)*B**(n-k)*d**k for k in range(n+1)),F(0))
    return P(5)/P(3)
def source_q(g,ends,times,oldq):
    assert len(ends)==len(times)==73 and times[24]=='1/64' and times[32]=='1/32' and times[-1]=='1'
    gl,gu=map(F,g);e=[tuple(map(F,x)) for x in ends];dl,du=e[24]
    assert 0<gl<=gu and all(0<lo<=hi for lo,hi in e) and max(gl,dl)<=min(gu,du)
    assert 0<e[32][0]<gu
    L=max(gl,dl);K=F(0);terms=[]
    for j in range(24,73):
        mass=e[j][1]-(e[j+1][0] if j<72 else F(0));assert mass>0
        w=weight(times[j],times[j+1] if j<72 else None);v=mass*w;K+=v;terms.append((mass,w,v))
    band=K/(16*L);split=F(256)*(F(7,10)-F(8,15)*e[32][0]/gu);q=min(F(oldq),band,F(512,3),split);assert 0<q<=F(oldq)
    return {'K':K,'lower64':L,'band':band,'split':split,'q':q,'terms':terms}
def corners(point,delta,gamma,q,lo,hi):
    p=list(map(F,point));delta=F(delta);g=[tuple(map(F,x)) for x in gamma];q=F(q);lo=list(map(F,lo));hi=list(map(F,hi))
    assert len(p)==5 and p[0]==1 and len(g)==len(lo)==len(hi)==4 and all(a<=b for a,b in g) and all(a<=b for a,b in zip(lo,hi)) and q>0
    D=lcm(*(x.denominator for x in [*p,delta,q,*lo,*hi,*[v for r in g for v in r]]));integer=lambda x:int(x*D)
    P=list(map(integer,p[1:]));A=integer(delta);Q=integer(q);G=[tuple(map(integer,r)) for r in g];L=list(map(integer,lo));U=list(map(integer,hi));values=[]
    for bits in product((0,1),repeat=4):
        X=[(U[j] if bits[j] else L[j])-P[j] for j in range(4)]
        N=A*D*D+D*sum(min(a*x,b*x) for (a,b),x in zip(G,X))-Q*sum(x*x for x in X);values.append((N,bits))
    n,bits=min(values);return F(n,D**3),bits
