#!/usr/bin/env python3
"""Exact Fraction certificates for the round31 special parameter choice."""
from fractions import Fraction as Q
n=Q(2); gamma=Q(1,2); beta=1-gamma; d=Q(1,10); tau=n*beta+d; delta=d/tau
assert delta==Q(1,11) and 0<delta<Q(1,3)
# B=2 sqrt(5)/(3*11^(1/10)); compare tenth powers exactly.
B10=Q(2,3)**10*Q(5)**5/Q(11)
assert B10==Q(3200000,649539) and B10>Q(7,6)**10
# Finite instance: a=s^18, s=1/10, so fractional powers are integral powers of s.
s=Q(1,10); a=s**18
u=s**20+Q(5)*s**18*(1-s**2)
F=s**11+Q(11,3)*s**9*(1-s**2)
assert u==Q(31,6250000000000000000) and F==Q(91,25000000000)
eta=Q(1,1000); t=u/(100*(1-eta)); assert t==Q(31,624375000000000000000)
# After cancelling the positive cube volume V: V delta^2 u=V(1-eta)t(1-delta)^2.
assert delta**2*u==(1-eta)*t*(1-delta)**2
# z=delta^(11/10)F has rational tenth power.  This implies z>(7/6)sqrt(t).
z10=delta**11*F**10; target10=Q(7,6)**10*t**5
assert 0<z10<1 and z10>target10
print({'B_to_10':str(B10),'B10_minus_(7_over_6)^10':str(B10-Q(7,6)**10),'a':str(a),'u':str(u),'F':str(F),'t':str(t),'z10_minus_target10':str(z10-target10),'status':'PASS'})
