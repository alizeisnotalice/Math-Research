"""Direct EPS128 rational-pole formula, integer outward arithmetic only.

Formula/rounding design is derived from accepted081; no numerical library or
production DD module is imported. No logarithm/constant/H is computed on import.
"""
from fractions import Fraction as F
from math import lcm,prod
from functools import lru_cache

PRECISION=768
S=1<<PRECISION
LOG_N=282
ATAN_N=216

def up(n,d):assert d>0;return -((-n)//d)

@lru_cache(maxsize=512)
def atanh(n,d):
    assert 0<=3*n<=d
    lo=S*n//d;hi=up(S*n,d);sqlo=lo*lo//S;sqhi=up(hi*hi,S);plo,phi=lo,hi;a=b=0
    for j in range(LOG_N):
        odd=2*j+1;a+=plo//odd;b+=up(phi,odd);plo=plo*sqlo//S;phi=up(phi*sqhi,S)
    odd=2*LOG_N+1
    return 2*a,2*b+up(9*S,4*odd*3**odd)

def logarithm(n,d=1):
    assert n>0 and d>0
    e=n.bit_length()-d.bit_length()
    if e>=0:
        if n<d<<e:e-=1
    elif n<<(-e)<d:e-=1
    a,b=(n,d<<e) if e>=0 else (n<<(-e),d)
    assert b<=a<2*b
    lo,hi=atanh(a-b,a+b);l2,u2=atanh(1,3)
    return (lo+e*l2,hi+e*u2) if e>=0 else (lo+e*u2,hi+e*l2)

def atan_inverse(q):
    a=b=0;power=q
    for j in range(ATAN_N):
        den=(2*j+1)*power;lo=S//den;hi=up(S,den)
        if j%2:a-=hi;b-=lo
        else:a+=lo;b+=hi
        power*=q*q
    tail=up(S,(2*ATAN_N+1)*power)
    if ATAN_N%2:a-=tail
    else:b+=tail
    return a,b

def pi_bounds():
    a,b=atan_inverse(5);c,d=atan_inverse(239)
    assert 476*239+480==480*239-476==114244
    lo,hi=16*a-4*d,16*b-4*c;assert 3*S<lo<hi<4*S
    return lo,hi

def add(acc,n,d,interval):
    if d<0:n,d=-n,-d
    assert d>0;a,b=interval
    if n<0:a,b=b,a
    acc[0]+=n*a//d;acc[1]+=up(n*b,d)

def values(point,fan):
    """3504 cones, direct pole weights and their four physical derivatives."""
    p=list(map(F,point));assert len(p)==5 and p[0]==1
    den=lcm(*(x.denominator for x in p));ints=[int(x*den) for x in p]
    assert fan['denominator']==144 and len(fan['vertices'])==280 and len(fan['simplices'])==3504
    k=[sum(a*b for a,b in zip(v,ints)) for v in fan['vertices']]
    assert len(set(k))==280 and 0 not in k
    # a=144*D*epsilon, epsilon fixed to1/128. This is direct EPS128,
    # without the EPS64(2x) scaling route used by frozen116.
    a=F(144*den,128);A,B=a.numerator,a.denominator
    sq=[A*A+B*B*x*x for x in k];logs=[logarithm(w) for w in sq]
    slopes=[[den*v for v in row[1:]] for row in fan['vertices']]
    sums=[[0,0] for _ in range(5)]
    for cone in fan['simplices']:
        ids=cone['vertices'];det=cone['det_numerator'];assert len(ids)==len(set(ids))==5 and det>0
        for j in ids:
            others=[m for m in ids if m!=j];diff=[k[j]-k[m] for m in others];D=prod(diff)
            add(sums[0],-det,2*D,logs[j])
            quotients=[D//v for v in diff]
            for axis in range(4):
                T=sum((slopes[j][axis]-slopes[m][axis])*q for m,q in zip(others,quotients))
                add(sums[axis+1],det*T,2*D*D,logs[j])
                add(sums[axis+1],-det*B*B*slopes[j][axis]*k[j],D*sq[j],(S,S))
    pil,piu=pi_bounds();num=den**4*S**4;dens=(4608*pil**5,4608*piu**5);out=[]
    for low,high in sums:
        endpoints=[F(num*v,d) for v in (low,high) for d in dens];out.append((min(endpoints),max(endpoints)))
    assert out[0][0]>0
    return out,dict(cones=3504,common_denominator=str(den),distinct_projections=280,minimum_integer_projection_gap=str(min(y-x for x,y in zip(sorted(k),sorted(k)[1:]))),epsilon='1/128')

def dyadic(ball,bits=512):
    scale=1<<bits;a,b=ball
    return [str(F((a*scale).__floor__(),scale)),str(F((b*scale).__ceil__(),scale))]
