#!/usr/bin/env python3
"""Independent rational enclosure of B_n(t) and (sqrt(A)+sqrt(B))^2."""
import sys
if hasattr(sys,"set_int_max_str_digits"): sys.set_int_max_str_digits(0)
from fractions import Fraction as Q
from math import comb, isqrt

def floorq(x, d): return Q(x.numerator*d//x.denominator,d)
def ceilq(x,d): return -floorq(-x,d)
def log_unit(x,N=250):
    z=(x-1)/(x+1)
    assert 0 <= z <= Q(1,3)
    s=2*sum((z**(2*j+1)/ (2*j+1) for j in range(N)),Q())
    return s,s+2*z**(2*N+1)/((2*N+1)*(1-z*z))
def log_bounds(t):
    k=0; x=t
    while x>=2: x/=2; k+=1
    l,u=log_unit(x); l2,u2=log_unit(Q(2))
    d=10**210
    return floorq(l+k*l2,d),ceilq(u+k*u2,d)
def sqrt_bounds(x):
    assert x>=0
    d=10**90
    m=isqrt(x.numerator*d*d//x.denominator)
    return Q(m,d),Q(m+1,d)
def certify(n,t):
    assert n>=1 and t>1
    q=1+1/(2*t); h=(1+1/t)/2
    beta=n-Q(n,2)*q**(n-1)
    terms=sum((Q(comb(n-1,j),j)*q**(n-1-j)*Q(-1,2)**j*(1-t**(-j)) for j in range(1,n)),Q())
    gamma=1-h**n-Q(n,2)*terms
    l,u=log_bounds(t)
    lo=gamma+beta*(l if beta>=0 else u)
    hi=gamma+beta*(u if beta>=0 else l)
    d=10**80; lo=floorq(lo,d); hi=ceilq(hi,d)
    assert 0<lo<=hi
    A=(1+1/t)**n
    sl,_=sqrt_bounds(A*lo); _,su=sqrt_bounds(A*hi)
    flo=floorq(A+lo+2*sl,d); fhi=ceilq(A+hi+2*su,d)
    return A,lo,hi,flo,fhi
def finite_certify(n,t,Bup):
    a=2**256; u=Q(1,a)
    assert a>2*t
    Aa=(1+1/t+u)**n
    logup=log_bounds(t)[1]
    delta=((1+u)**n-1)*(1+n*logup)+u**n
    Eup=Bup+delta
    _,tau=sqrt_bounds(Aa*Eup)
    finite=(1+Aa/tau)*(tau+Eup)
    return Aa,delta,Eup,tau,finite

if __name__=='__main__':
    import json, pathlib
    data=json.loads((pathlib.Path(__file__).parent/'table.json').read_text())
    for row in data:
        n=row['n']; t=Q(row['t']); A,lo,hi,fl,fu=certify(n,t)
        assert (str(lo),str(hi),str(fl),str(fu)) == tuple(row[k] for k in ('B_lower','B_upper','F_lower','F_upper'))
        U=Q(row['U']); assert fu<U
        Aa,delta,Eup,tau,finite=finite_certify(n,t,hi)
        assert tuple(map(str,(Aa,delta,Eup,tau,finite))) == tuple(row[k] for k in ('finite_A','finite_delta_upper','finite_E_upper','finite_tau','finite_bound'))
        assert finite<U
        if n==7:
            assert t==Q(3391749,524288) and U==Q(4910992105193,200000000000)
    print('PASS:',len(data),'independent exact rational certificates')
