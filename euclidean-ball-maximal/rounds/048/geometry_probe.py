#!/usr/bin/env python3
"""Two-source geometry: exact n1 areas, Decimal envelope checks, exact no-go box."""
from fractions import Fraction as F
from decimal import Decimal as D,localcontext,ROUND_FLOOR
from itertools import product
import json

def area_outside(xa,xb,za,zb,R):
 if xa>=xb or za>=zb:return F(0)
 square=lambda v:max(F(0),v)**2
 strip=sum((sign*(square(zb+shift)-square(za+shift))/2 for shift,sign in ((R-xa,1),(R-xb,-1),(-R-xa,-1),(-R-xb,1))),F(0))
 result=(xb-xa)*(zb-za)-strip
 nodes=sorted({za,zb}|{q for q in (xa-R,xb-R,xa+R,xb+R) if za<q<zb})
 def length(z):return max(F(0),min(xb,z-R)-xa)+max(F(0),xb-max(xa,z+R))
 other=sum(((v-u)*length((u+v)/2) for u,v in zip(nodes,nodes[1:])),F(0))
 assert result==other and 0<=result<=(xb-xa)*(zb-za)
 return result

def n1_checks():
 out=[]
 for batch,(L,den) in enumerate(((3,8),(8,32),(24,128))):
  maxima=[];checks=0
  for l in range(1,L+1):
   r=F(1,2**l);best=(F(0),F(0))
   for i in range(2*den+1):
    d=F(i,den);xa=max(F(-1),d-1);xb=min(F(1),d+1);za=-r;zb=min(r,d-1)
    T=area_outside(xa,xb,za,zb,F(1))/(4*r)
    assert 0<=T<=r/2 # sharper n1 translation-loss constant beta_1=1/2
    if d<=1-r or d>=2:assert T==0
    if T>best[0]:best=(T,d)
    checks+=1
   maxima.append(dict(l=l,Tmax_grid=str(best[0]),d=str(best[1])))
  out.append(dict(batch=batch,max_gap=L,source_distance_denominator=den,exact_area_checks=checks,maximum_by_gap=maxima))
 return out

def envelope(n,prec):
 with localcontext() as ctx:
  ctx.prec=prec;nn=D(n);k=D(2).ln()/nn;q=(-k).exp();root=nn.sqrt();c=D(3)/64;A=(-c*nn).exp()
  prefix=D(0)
  for l in range(1,n+1):
   t=(-k*l).exp();prefix+=min((-(nn*(1-t*t))/16).exp(),root*t)
  cross=((root/A).ln()/k).to_integral_value(rounding=ROUND_FLOOR);m=int(cross);start=max(n,m)
  tail=A*max(0,m-n)+root*(-k*(start+1)).exp()/(1-q)
  integral=A*nn/D(2).ln()*(1+D('0.5')*nn.ln()+c*nn)
  assert tail<=integral+D('1e-40') and prefix+tail<83
  return prefix+tail,prefix,tail

def high_dimension_checks():
 out=[];scalar=0
 batches=((1,2,4),(16,32,64),(128,1024,4096))
 for batch,ns in enumerate(batches):
  for n in ns:
   val,prefix,tail=envelope(n,80);low,*_=envelope(n,50);assert abs(val-low)<D('1e-35')
   with localcontext() as ctx:
    ctx.prec=80
    for l in sorted({1,2,max(1,n//2),n,2*n}):
     t=(-D(2).ln()*l/n).exp();s=1-t*t;bound=(-D(n)*s/16).exp()
     for i in range(65):
      d=D(i)/32
      if not d or d>=2:surrogate=D(0)
      else:
       inter=(D(n)/2*(1-d*d/4).ln()).exp()
       if d*d>=s:cap=D(1)
       else:
        a=(s-d*d)/(2*d)
        cap=D(0) if a>=t else (D(n)/2*(1-a*a/(t*t)).ln()).exp()
       surrogate=inter*cap
      assert surrogate<=bound+D('1e-60');scalar+=1
   out.append(dict(batch=batch,n=n,envelope_upper_numeric=str(val),prefix_numeric=str(prefix),tail_numeric=str(tail),precision_digits=[50,80],cross_precision_error_less_than='1e-35'))
 return out,scalar

def uniform_mass(x,r):
 return sum((height*max(F(0),min(b,x+r)-max(a,x-r)) for a,b,height in ((F(-1,12),F(1,12),F(3,2)),(F(10),F(11),F(3,4)))),F(0))
def uniform_max(x):
 endpoints=(F(-1,12),F(1,12),F(10),F(11));radii={abs(x-a) for a in endpoints}-{F(0)}
 return max(uniform_mass(x,r)/(2*r) for r in radii)

def relaxed_box_checks():
 rows=[]
 for batch,js in enumerate(((8,9,12),(16,32,64),(128,256,512))):
  for j in js:
   R=F(2)**(2-j);r=R/2
   for y,s,d,t in product((F(-1,48),F(1,48)),(F(3,4),F(7,8)),(F(3,4),F(7,8)),(F(-1,2),F(-3,8))):
    x=y+s*R;v=y+d*R;z=y+t*R
    assert abs(x)<F(1,12) and abs(z)<F(1,12) and abs(y)<F(1,12) and abs(v)<F(1,12)
    assert abs(x-y)<R and abs(x-v)<R and abs(z-y)<=r and abs(z-v)>R and abs(x-z)>R
    for q in (x,z):
     assert uniform_max(q)==F(3,2)
     trace=[uniform_mass(q,F(2)**(2-k))/(2*F(2)**(2-k)) for k in range(1,6)]
     assert trace[:4]==[F(1,16),F(1,8),F(1,4),F(1,2)] and trace[4]>F(1,2)
   volume=F(1,24)*(R/8)**3
   charge=4*F(9,4)*volume/((2*R)**2*(2*r))
   assert charge==F(3,16384)
   rows.append(dict(batch=batch,j=j,k=j+1,per_scale_lower=str(charge),jacobian_volume=str(volume)))
 # Scalar constants in the all-n proof are certified rationally.
 c=F(3,64);elow=F(8,3);log2low=F(2,3)
 assert 1/(c*elow)==8 and 4/(c*elow**2)==12 and 12**3<42**2
 assert F(64,3)+(8+21+12)/log2low==F(497,6)<83
 return rows

def main():
 one=n1_checks();high,scalar=high_dimension_checks();box=relaxed_box_checks()
 result=dict(status='passed',n1=one,high_dimension=high,scalar_cap_checks=scalar,relaxed_box=box,
  analytic_all_dimension_sum_bound=83,analytic_rational_upper='497/6',
  high_dimension_scope='50/80-digit Decimal evaluation of proved upper envelopes, not interval volume computation',
  exact_n1_area_checks=sum(r['exact_area_checks'] for r in one),
  relaxed_input='(3/2)1_(-1/12,1/12)+(3/4)1_(10,11)',actual_outer_cost='0',relaxed_outer_lower='3*(D-8)/16384 for D>=9',
  main_theorem_proved=False)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
