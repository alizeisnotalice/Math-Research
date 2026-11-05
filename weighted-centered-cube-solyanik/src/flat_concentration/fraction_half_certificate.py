#!/usr/bin/env python3
"""Exact beta=1/2, u=1/8 certificate in common denominator 104."""
from fractions import Fraction as Q
p,q,r,N=Q(2),Q(3),Q(10),Q(104)
a=Q(1,26);u=Q(1,8);beta=Q(1,2)
# p,q,r are the positive square-root numerators for a,u-a,1-a.
assert p*p/N==a and q*q/N==u-a and r*r/N==1-a
assert 0<a<u/2 and 0<u<1
omega=(p+q)/(p+r)
# K after cancellation of the common sqrt(N) factors:
# beta*((1/p-1/q)*(p+r)-(p+q)*(1/p-1/r)).
K=beta*((Q(1,p)-Q(1,q))*(p+r)-(p+q)*(Q(1,p)-Q(1,r)))
assert omega==Q(5,12) and K==0
print({'u':str(u),'a':str(a),'omega':str(omega),'K':str(K),'status':'PASS'})
