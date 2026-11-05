#!/usr/bin/env python3
"""Exact rational checks for the Round32 Cartesian-grid lower construction."""
from fractions import Fraction as Q

# Lambda for gamma=1/5, beta=4/5: compare fifth powers exactly.
gamma=Q(1,5); beta=1-gamma
lambda5=Q(2)**4 * Q(4,5)**4 * Q(1,5)
assert lambda5==Q(4096,3125) and lambda5>1

# Finite n=2 instance.  h includes the strict missing factor 1-eta.
n=Q(2); gamma=Q(1,20); beta=1-gamma; k=Q(2); ell=Q(3,38); eta=Q(1,76)
h=(1-eta)*(Q(2)/ell)**(n-1); T=k*(1+ell)+1
assert h==25 and T==Q(60,19) and h*ell/Q(2)==1-eta
band_fraction=k*ell/T
assert band_fraction==Q(1,20)
# J=(k ell/T)^(n-1) h^beta, so J^20 is rational.
J20=band_fraction**20 * h**19
assert J20==Q(5**18,2**40) and J20>2
# Hence J>2^(1/20); this is the positive leading coefficient condition.
print({'Lambda_1_5_to_5':str(lambda5),'h':str(h),'T':str(T),'band_fraction':str(band_fraction),'J_to_20':str(J20),'J20_minus_2':str(J20-2),'status':'PASS'})
