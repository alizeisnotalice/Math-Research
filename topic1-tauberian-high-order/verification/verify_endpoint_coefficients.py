#!/usr/bin/env python3
"""Exact polynomial certificate for all 0 < epsilon < 1, with no symbolic library."""
from fractions import Fraction as F
import json,time

indices=(0,1,3,4,5,6,8,9)
def add(a,b): return [a[i]+b[i] for i in range(4)]
def scale(a,c): return [x*c for x in a]
def mul(a,b):
 out=[F(0)]*4
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<4:out[i+j]+=x*y
 return out

def main():
 start=time.time(); counts={'positive':0,'negative':0,'identically_zero':0,'endpoint_root_zero':0,'endpoint_root_one':0,'terms':0}; roots=[]
 def positive_cube(z):
  # z is affine, and an affine function on (0,1) has the sign of its midpoint
  # if endpoint signs do not strictly disagree.
  a,b=z; counts['terms']+=1
  if b:
   root=-a/b
   assert not (0<root<1), ('positive-part sign changes',z,root)
   if root==0:counts['endpoint_root_zero']+=1
   if root==1:counts['endpoint_root_one']+=1
  mid=a+b/2
  if mid>0:
   counts['positive']+=1; v=[a,b,F(0),F(0)]; return mul(mul(v,v),v)
  if mid<0:counts['negative']+=1;return [F(0)]*4
  assert a==b==0
  counts['identically_zero']+=1;return [F(0)]*4
 def difference(x,y,k):return [x[0]-y[0]-3*k,x[1]-y[1]]
 def antiderivative(z,r):
  x=positive_cube([z[0]+r,z[1]])
  x=add(x,scale(positive_cube(z),-2))
  x=add(x,positive_cube([z[0]-r,z[1]]))
  return scale(x,F(1,6)/r)
 polys={F(1):[F(0)]*4,F(1,4):[F(0)]*4}
 for j in indices:
  r=F(1) if j in (0,1,8,9) else F(1,4)
  a=[F(j,4),F(0)];b=[F(j,4),F(1,4)]
  for k in (-1,0,1):
   for h in indices:
    c=[F(h,4),F(0)];d=[F(h,4),F(1,4)]
    # Integral of T(x-y) over x in [a,b], y in [c+3k,d+3k].
    piece=[F(0)]*4
    for p,q,sign in ((b,c,1),(a,c,-1),(b,d,-1),(a,d,1)):
     piece=add(piece,scale(antiderivative(difference(p,q,k),r),sign))
    polys[r]=add(polys[r],piece)
 assert polys[F(1)]==[F(0),F(0),F(19,32),F(0)]
 assert polys[F(1,4)]==[F(0),F(0),F(1,4),F(-1,48)]
 out={'interval':'0 < epsilon < 1','source_and_target_intervals':len(indices),'period_shifts':[-1,0,1],'positive_part_sign_checks':counts,'all_affine_positive_part_roots_outside_open_interval':True,'coefficients_order':'constant, epsilon, epsilon^2, epsilon^3','energy_by_radius':{str(r):list(map(str,p)) for r,p in polys.items()},'all_expected_coefficients_match':True,'elapsed_seconds':time.time()-start}
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
