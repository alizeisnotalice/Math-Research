"""Outward rational scalar arithmetic only; no computation on import."""
from fractions import Fraction as F

BITS=384
LOG_BITS=512
LOG_N=192

def interval(x):
    if isinstance(x,(tuple,list)):return tuple(map(F,x))
    q=F(x);return q,q

def add(a,b):
    a,b=interval(a),interval(b);return a[0]+b[0],a[1]+b[1]

def sub(a,b):
    a,b=interval(a),interval(b);return a[0]-b[1],a[1]-b[0]

def mul(a,b):
    a,b=interval(a),interval(b);v=[x*y for x in a for y in b];return min(v),max(v)

def div(a,b):
    a,b=interval(a),interval(b);assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))

def power(a,n):
    a=interval(a);assert a[0]>=0 and type(n) is int and n>0
    return a[0]**n,a[1]**n

def ceildiv(n,d):assert d>0;return -((-n)//d)

def root_scalar(q,n):
    q=F(q);assert q>=0 and n>=2
    scale=1<<BITS;v=q.numerator*scale**n//q.denominator
    lo,hi=0,1<<((v.bit_length()+n-1)//n)
    while hi-lo>1:
        m=(lo+hi)//2
        if m**n<=v:lo=m
        else:hi=m
    low,high=F(lo,scale),F(lo+1,scale)
    assert low**n<=q<high**n
    return low,high

def root(a,n):
    a=interval(a);assert a[0]>=0
    return root_scalar(a[0],n)[0],root_scalar(a[1],n)[1]

def _atanh(n,d):
    assert 0<=3*n<=d
    S=1<<LOG_BITS;lo,hi=S*n//d,ceildiv(S*n,d);sl,sh=lo*lo//S,ceildiv(hi*hi,S);pl,ph=lo,hi;a=b=0
    for j in range(LOG_N):
        odd=2*j+1;a+=pl//odd;b+=ceildiv(ph,odd);pl=pl*sl//S;ph=ceildiv(ph*sh,S)
    odd=2*LOG_N+1;tail=ceildiv(9*S,4*odd*3**odd)
    return 2*a,2*b+tail

def log_scalar(q):
    q=F(q);assert q>0;n,d=q.numerator,q.denominator;e=n.bit_length()-d.bit_length()
    if e>=0:
        if n<d<<e:e-=1
    elif n<<(-e)<d:e-=1
    a,b=(n,d<<e) if e>=0 else (n<<(-e),d);assert b<=a<2*b
    lo,hi=_atanh(a-b,a+b);l2,u2=_atanh(1,3)
    lo,hi=(lo+e*l2,hi+e*u2) if e>=0 else (lo+e*u2,hi+e*l2)
    S=1<<LOG_BITS;return F(lo,S),F(hi,S)

def log(a):
    a=interval(a);assert a[0]>0
    return log_scalar(a[0])[0],log_scalar(a[1])[1]

def pack(a):
    a=interval(a);assert a[0]<=a[1]
    return dict(lower=str(a[0]),upper=str(a[1]),lower_display=float(a[0]),upper_display=float(a[1]))
