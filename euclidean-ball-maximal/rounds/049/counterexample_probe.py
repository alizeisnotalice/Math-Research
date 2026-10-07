#!/usr/bin/env python3
"""Exact mixed-density global-comparison obstruction and activation budgets."""
from fractions import Fraction as F
import json

def mass(x,r,m):
 return F(1,2)*max(F(0),min(F(1,2),x+r)-max(F(0),x-r))+(m if abs(x)<r else 0)+(F(3,4)-m if abs(x-10)<r else 0)
def maximal(x,m):
 if x in (0,10):return None
 radii={abs(x),abs(x-F(1,2)),abs(x-10)}-{F(0)}
 def critical(r):
  return (F(1,2)*max(F(0),min(F(1,2),x+r)-max(F(0),x-r))+(m if abs(x)<=r else 0)+(F(3,4)-m if abs(x-10)<=r else 0))/(2*r)
 return max([F(1,2) if 0<x<F(1,2) else F(0)]+[critical(r) for r in radii])
def R(l):return F(2)**(2-l)
def trace(x,m,D):return [mass(x,R(l),m)/(2*R(l)) for l in range(1,D+1)]
def first(t,b):return next((j for j,g in enumerate(t,1) if g>b),None)
def main():
 rows=[];checks=0;pack=0;smooth=0
 for batch,ks in enumerate(((8,10,12),(16,24,32),(48,64,96))):
  for k in ks:
   m=R(k)/5;D=k+2;X=F(11,12)*m
   for z in (m/2,5*m/8,3*m/4):
    tr=trace(z,m,D);assert 1<maximal(z,m)<2 and first(tr,F(1,8))>=2 and first(tr,F(1,2))<=D
    for beta in (F(7,20),F(9,25),F(37,100)):
     assert first(tr,beta)==k
     c=4*beta-1
     for l in range(5,k):
      a=c*R(l)-2*m;b=2*c*R(l)-2*m
      assert z-R(l)<a<b<=z+R(l)
      assert (b-a)/(2*R(l))==c/2>=F(1,5)
      for t in (F(1,8),F(1,2),F(7,8)):
       q=a+t*(b-a);tq=trace(q,m,D)
       assert first(tq,beta)==l
       checks+=1
   # Entire original central E has no coarse exits j<k for beta>=7/20.
   for z in [-m/2,-3*m/8,-m/4,m/3,m/2,3*m/4,m]:
    tr=trace(z,m,D)
    assert max(tr[:k-1])<F(7,20)
   # Exact central maximal formula and far band J1, sampled endpoints/interiors.
   for z in [-2*m,-m,-m/2,-3*m/8,-m/4,m/4,m/3,m/2,3*m/4,m,2*m]:
    expected=m/(2*abs(z)) if z<0 else F(1,2)+m/(2*z)
    assert maximal(z,m)==expected
   w=F(3,4)-m
   for z in [10-w/2,10-3*w/8,10-w/4,10+w/4,10+3*w/8,10+w/2]:
    assert maximal(z,m)==w/(2*abs(z-10)) and trace(z,m,D)[0]>F(1,8)
   # Integral over z core of g_k is 61m/640; averaged c/2 over beta is 11/625.
   zint=F(7,20)*(m/4)+((3*m/4)**2-(m/2)**2)/(40*m)
   beta_lo=F(7,20);beta_hi=F(37,100)
   beta_factor=4*((beta_hi**2-beta_hi/2)-(beta_lo**2-beta_lo/2))
   bound=(k-5)*zint*beta_factor/X
   assert zint==F(61,640)*m and beta_factor==F(11,625)
   assert bound==F(183,100000)*(k-5)
   # For each actual single exit j, every finite tail packing is explicit.
   for j in (2,k,D):
    for L in (j,D,D+8):
     weights=[F(2)**(j-l) if l>=j else F(0) for l in range(1,L+1)]
     assert sum(weights)==2-F(2)**(j-L) and max(weights)<=1
     assert sum(2*w for w in weights)<=4
     pack+=1
   h=m/1000
   assert 2*(m+h)<3*m
   for z in (m/2,3*m/4):
    lower=F(1,4)+(z-h+2*m)/(4*(z+h))
    assert lower>F(29,25) and maximal(z-h,m)<F(8,5)
    for tau in (-h,h):
     tz=trace(z-tau,m,D)
     assert first(tz,F(7,20))==first(tz,F(37,100))==k
     assert first(tz,F(1,8))>=2 and first(tz,F(1,2))<=D
   for beta in (F(7,20),F(37,100)):
    c=4*beta-1
    for l in range(5,k):
     a=(c+F(1,50))*R(l)-2*m;b=(2*c-F(1,50))*R(l)-2*m
     assert (b-a)/(2*R(l))>=F(9,50)>F(1,6)
     for q in (a,b):
      for tau in (-h,h):
       tq=trace(q-tau,m,D)
       assert first(tq,beta)==l
       assert tq[l-1]==F(1,4)+(q-tau+2*m)/(4*R(l))
       smooth+=1
   assert F(7,20)*(m/4)*F(1,6)*(F(2,25))/(3*m)==F(7,18000)
   rows.append(dict(batch=batch,k=k,D=D,m=str(m),X=str(X),beta_interval=['7/20','37/100'],coarse_layers=k-5,mean_Bglobal_over_X_lower=str(bound),original_coarse_activation_on_test_block='zero'))
 print(json.dumps(dict(status='passed',cases=rows,batch_counts=[3,3,3],exact_global_label_checks=checks,activation_scalar_checks=pack,smooth_extreme_shift_checks=smooth,smooth_Bglobal_ratio_lower="7*(k-5)/18000",analytic_mean_Bglobal_over_X_lower='183*(k-5)/100000',scope='Finite rational checks supplement explicit all-k proof; do not measure actual compatible cost.',main_theorem_proved=False),indent=2))
if __name__=='__main__':main()
