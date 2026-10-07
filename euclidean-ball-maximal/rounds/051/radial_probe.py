#!/usr/bin/env python3
"""Actual two-dimensional clocks and radial-label loss, using scaled rational coordinates."""
from fractions import Fraction as F
import json
A=[((F(0),F(0)),F(17,400)),((F(6,5),F(0)),F(11,200)),((F(67,100),F(0)),F(1,8000)),((F(100),F(0)),F(7219,8000))]
def d2(a,b):return sum((x-y)**2 for x,y in zip(a,b))
def r2(i):return F(2)**(5-i)
def trace(x):return [4*sum((w for a,w in A if d2(x,a)<r2(i)),F(0))/r2(i) for i in range(1,8)]
def first(t,b):return next((i for i,v in enumerate(t,1) if v>b),None)
def maximal(x):
 ds=sorted({d2(x,a) for a,w in A});assert ds[0]>0
 return max(4*sum((w for a,w in A if d2(x,a)<=r),F(0))/r for r in ds)
def main():
 x=(F(13,20),F(0));z=(-F(2,5),F(0));tx=trace(x);tz=trace(z)
 assert sum(w for a,w in A)==1
 assert maximal(x)==F(5,4) and maximal(z)==F(17,16)
 assert first(tx,F(1,8))==4 and first(tz,F(1,8))==5
 assert first(tx,F(1,2))==6 and first(tz,F(1,2))==7
 # sqrt(2) is bracketed rationally: this is a certificate, not a floating test.
 lo=F(1414213,1000000);hi=F(1414214,1000000)
 assert lo**2<2<hi**2
 cases=[];checks=0
 for batch,eps in enumerate((F(1,10000),F(1,100000),F(1,1000000))):
  for sx in (-1,1):
   for sy in (-1,1):
    xx=(x[0]+sx*eps,sy*eps)
    for sz in (-1,1):
     for sw in (-1,1):
      zz=(z[0]+sz*eps,sw*eps)
      assert 1<maximal(xx)<2 and 1<maximal(zz)<2
      assert d2(xx,zz)>1 and d2(xx,(0,0))<F(1,2) and d2(zz,(0,0))<F(1,4)
      for beta in (F(3,8),F(19,50)):
       assert first(trace(xx),beta)==5 and first(trace(zz),beta)==7
       checks+=1
      assert first(trace(xx),F(1,2))==6 and first(trace(zz),F(1,2))==7
      # A rectangle enclosing u=sqrt(2)*xx is within distance sqrt(.11) of v.
      maxdx=max(abs(F(6,5)-lo*xx[0]),abs(F(6,5)-hi*xx[0]));maxdy=hi*abs(xx[1])
      assert maxdx**2+maxdy**2<F(11,100)
  cases.append(dict(batch=batch,box_halfwidth=str(eps),beta_over_alpha=['3/8','19/50'],x_label=5,z_label=7,radial_image_MP_over_alpha_lower=2,strict=True))
 print(json.dumps(dict(status='passed',n=2,alpha='8/pi',coordinate_unit='R5=1/sqrt(32)',D=7,atoms=[dict(position=[str(v) for v in a],mass=str(w)) for a,w in A],center_MP_over_alpha=['5/4','17/16'],J=[4,5],half_threshold_labels=[6,7],clock_checks=checks,cases=cases,scope='Actual first-exit counterexample to radial monotonicity. Finite exact tests support the strict-neighborhood proof; no divergent weak-type conclusion.'),indent=2))
if __name__=='__main__':main()
