#!/usr/bin/env python3
"""Exact all-clock and full-radius tests for one-dimensional backward shear."""
from fractions import Fraction as F
import json
from decimal import Decimal as Dec,localcontext

def R(i):return F(2)**(2-i)
def atoms(j):
 s=F(2)**(5-j)
 return [(F(-1,2)+F(3,100)*s,F(19,100)),(-F(3,200)*s,F(1,25)*s),(F(21,200)*s,F(3,50)*s),(F(1,2)+F(1,50)*s,F(19,100)),(F(10),F(31,50)-s/10)]
def trace(A,x,D):return [sum((w for a,w in A if abs(x-a)<R(i)),F(0))/(2*R(i)) for i in range(1,D+1)]
def first(t,b):return next((i for i,g in enumerate(t,1) if g>b),None)
def maximal(A,x):
 ds=sorted({abs(x-a) for a,w in A});assert ds[0]>0
 return max(sum((w for a,w in A if abs(x-a)<=r),F(0))/(2*r) for r in ds)
def main():
 cases=[];checks=0
 for batch,js in enumerate(((5,6,8),(12,20,32),(48,80,128))):
  for j in js:
   s=F(2)**(5-j);A=atoms(j);D=j+2;k=j+1;y=F(21,200)*s;eps=s/10000;h=eps/10
   assert sum(w for a,w in A)==1 and min(w for a,w in A)>0
   assert maximal(A,0)==F(4,3) and maximal(A,F(13,100)*s)==F(6,5)
   assert first(trace(A,F(1,40)*s,D),F(3,8))==3
   for x in (-eps,eps):
    for z in (F(13,100)*s-eps,F(13,100)*s+eps):
     assert abs(x-z)>R(j) and abs(x-y)<R(j) and abs(z-y)<R(k)
     for v in (x,z):
      mv=maximal(A,v);dmin=min(abs(v-a) for a,w in A)
      assert F(11,10)<mv<F(3,2)
      assert mv*dmin/(dmin+h)>1 and mv*dmin/(dmin-h)<2
     # Arbitrary smoothing of every source by shifts in [-h,h].
     for shift in (-h,h):
      shifted=[(a+shift,w) for a,w in A]
      tx=trace(shifted,x,D);tz=trace(shifted,z,D)
      assert first(tx,F(1,8))>=2 and first(tz,F(1,8))>=2
      assert first(tx,F(1,2))==first(tz,F(1,2))==D
      for beta in (F(3,8),F(379,1000)):
       assert first(tx,beta)==j and first(tz,beta)==k
       # Shear source shifts need not equal any particular background shift.
       for source_shift in (-h,h):
        q=x+z-y-source_shift;tq=trace(shifted,q,D)
        assert first(tq,beta)==3
        checks+=1
   fee=F(4,1)*(eps**2)*(F(3,50)*s)/(2*R(j)*2*R(k))*((F(379,1000)-F(3,8))/F(1,4))
   assert fee==F(12,9765625000)*s
   cases.append(dict(batch=batch,j=j,k=k,D=D,s=str(s),shear_at_center=str(s/40),backward_jump=j-3,box_halfwidth=str(eps),smoothing_radius=str(h),beta_interval=['3/8','379/1000'],positive_uncovered_outer_fee_lower=str(fee)))
 activation_checks=0
 for beta in (F(1,4),F(9,32),F(5,16),F(11,32),F(3,8),F(13,32),F(7,16),F(15,32),F(1,2)):
  for depth in (5,16,64):
   bound=sum((min(2*beta,2*F(2)**(-r)) for r in range(depth)),F(0))
   assert bound<=4*beta+1<=3
   activation_checks+=1
 assert 4*((2*F(1,2)**2+F(1,2))-(2*F(1,4)**2+F(1,4)))==F(5,2)
 factors=[]
 for batch,ns in enumerate(((1,2,3),(8,32,64),(128,1024,4096))):
  for n in ns:
   vals=[]
   for prec in (50,80):
    with localcontext() as ctx:
     ctx.prec=prec;t=(-Dec(2).ln()/n).exp();v=((1+t)/2)**n
     assert v<=Dec('0.75')+Dec('1e-45');vals.append(v)
   assert abs(vals[0]-vals[1])<Dec('1e-40')
   factors.append(dict(batch=batch,n=n,normalized_support_factor=str(vals[1]),precision_digits=[50,80]))
 print(json.dumps(dict(activation_scalar_checks=activation_checks,deep_factor_checks=factors,deep_factor_scope='High-precision scalar checks; all-n bound follows from concavity, not numerics.',status='passed' ,batch_counts=[3,3,3],cases=cases,strict_clock_extreme_checks=checks,uncovered_fee_lower='12*s/9765625000',scope='Strict actual-edge counterexample, not divergence of full O/X.',main_theorem_proved=False),indent=2))
if __name__=='__main__':main()
