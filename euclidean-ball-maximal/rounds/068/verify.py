#!/usr/bin/env python3
"""Exact certificates for the positive-kernel exterior obstruction."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys

def log_inverse_interval(x,terms=96):
 # -log(x) = 2 atanh((1-x)/(1+x)), with positive geometric tail.
 t=(1-x)/(1+x)
 lo=2*sum((t**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
 tail=2*t**(2*terms+1)/F(2*terms+1)/(1-t*t)
 return lo,lo+tail

def radial_intervals(n,z,terms=96):
 a=n//2
 # Direct integral of a*t^(a-1)/(1-t), expanded as a positive series.
 lo=a*sum((z**(a+k)/F(a+k) for k in range(terms)),F(0))
 hi=lo+a*z**(a+terms)/F(a+terms)/(1-z)
 # Independent polynomial division plus a rigorously bracketed logarithm.
 l,h=log_inverse_interval(1-z,terms)
 poly=sum((z**k/F(k) for k in range(1,a)),F(0))
 other_lo=a*(l-poly);other_hi=a*(h-poly)
 assert max(lo,other_lo)<=min(hi,other_hi)
 assert hi-lo<F(1,10**9) and other_hi-other_lo<F(1,10**9)
 return lo,hi

def floor40(x):
 scale=2**40
 return F(x.numerator*scale//x.denominator,scale)

def main():
 batches=[]
 for batch,ms in enumerate(((7,8,9),(10,11,12),(13,14,16))):
  rows=[]
  for m in ms:
   n=2**m;U=F(8*m,n);t=1-F(2*m,n)
   assert 0<U<1 and 0<t<1
   # t < (U/2)^(2/n), so the exact lower integral exceeds 2.
   assert t**(n//2)<U/2
   assert (U/2)/(1-t)==2
   eta=F(1,10*n);inner=(1-eta)**n;outer=(1+eta)**n
   assert inner>=F(9,10) and outer<=F(10,9)
   annulus=F(9,10)/U-F(10,9)
   assert annulus>0
   rows.append(dict(n=n,m=m,U=str(U),smooth_exterior_level_volume_lower_over_unit_ball=str(annulus),actual_inner_volume_lower=str(floor40(inner/U)),exterior_energy_lower_after_exception_4=str(4*max(F(0),annulus-4))))
  radial=0
  for n in ((2,4,8),(10,16,24),(32,48,64))[batch]:
   for z in (F(1,4),F(1,2),F(3,4)):
    lo,hi=radial_intervals(n,z)
    assert lo>0;radial+=1
  # Rational (deliberately weaker) n^2 coefficients for arbitrary exceptions.
  exceptions=[]
  for C in (1,4,16):
   K=4*(C+2);volume=F(9,10)*K-F(10,9)-C
   assert volume>0
   coeff=volume/F(64*K**4)
   exceptions.append(dict(exception_budget=C,K=K,energy_lower_coefficient_of_n_squared=str(coeff)))
  batches.append(dict(batch=batch,dyadic_dimensions=rows,independent_radial_integrals=radial,exception_certificates=exceptions))
 path=Path(__file__)
 print(json.dumps(dict(round=68,status='passed',python=sys.version.split()[0],batches=batches,sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest()},main_theorem_proved=False,scope='Positive potential on the source exterior, not integration restricted to actual maximal tasks.'),indent=2))
if __name__=='__main__':main()
