#!/usr/bin/env python3
"""Whole-ball radial-depth experiments and exact rational selection obstruction.
No numerical receipt is an interval-arithmetic proof. Analytic conclusions and
scope are stated separately in the report.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
from math import factorial,comb,cos,pi
import json,time,hashlib,argparse

def dec(x):return D(x.numerator)/D(x.denominator) if isinstance(x,F) else D(x)

def gauss(N):
 nodes=[];tol=D(10)**(-95)
 for i in range(1,N+1):
  x=D(str(cos(pi*(i-.25)/(N+.5))))
  for iteration in range(12):
   p0=D(1);p1=x
   for k in range(2,N+1):p0,p1=p1,((2*k-1)*x*p1-(k-1)*p0)/k
   derivative=N*(x*p1-p0)/(x*x-1);step=p1/derivative;x-=step
   if abs(step)<tol:break
  assert abs(step)<tol
  nodes.append((x,D(2)/((1-x*x)*derivative*derivative)))
 return nodes

def whole_ball(n,N,nodes):
 r=D('0.5');a=1-r;m=(n-1)//2
 normal=F(factorial(2*m+1),2**(2*m+1)*factorial(m)**2)
 coefficients=[dec(normal*F((-1)**k*comb(m,k)*2**(m-k),m+k+1)) for k in range(m+1)]
 def cdf(t):
  if t<=-1:return D(0)
  if t>=1:return D(1)
  if t>0:return 1-cdf(-t)
  u=1+t;p=D(0)
  for c in reversed(coefficients):p=p*u+c
  return u**(m+1)*p
 def outside_probability(rho):
  if rho<=a:return D(0)
  s=(r*r+rho*rho-1)/(2*rho)
  result=cdf(s/r)-cdf(s-rho)/(r**n)
  assert -D('1e-85')<=result<=1+D('1e-85')
  return max(D(0),min(D(1),result))
 def integrate(b,depth=False):
  if b<=a:return D(0)
  center=(a+b)/2;half=(b-a)/2;total=D(0)
  for t,w in nodes:
   rho=center+half*t;f=n*rho**(n-1)*outside_probability(rho)
   if depth:f*=(-n*rho.ln()/D(2).ln())
   total+=w*f
  return half*total
 outer=integrate(D(1));mean=integrate(D(1),True)/outer
 thresholds=(1,2,4,8,16);tails={str(t):integrate((-D(t)*D(2).ln()/n).exp())/outer for t in thresholds}
 partial_count=n;repeat=sum((integrate((-D(t)*D(2).ln()/n).exp())/outer for t in range(partial_count)),D(0))
 remainder=F(0) if partial_count==n else F(1,2**(partial_count-1))
 assert mean<1/D(2).ln()+D('1e-35')
 assert all(value<=D(2)**(-int(t))+D('1e-35') for t,value in tails.items())
 assert repeat<=2+D('1e-35')
 if n==1:
  assert abs(outer-D(1)/8)<D('1e-80')
  exactmean=(D(1)-D(1)/(D(2).ln()*2)) # r=1/2 exact integral
  assert abs(mean-exactmean)<D('1e-75')
 return dict(n=n,r_over_R='1/2',canonical_gap=n,quadrature_nodes=N,outer_probability=str(outer),conditional_depth_mean=str(mean),conditional_depth_tail={t:str(v) for t,v in tails.items()},conditional_repetition_partial=str(repeat),repetition_partial_terms=partial_count,analytic_uncomputed_repetition_tail_bound=str(remainder))

def exact_rays():
 batches=((1,2,4),(8,16,32),(64,128,256));rows=[];checks=0
 sources=[(F(-7,16),F(0)),(F(-1,4),F(0)),(F(1,4),F(0)),(F(-3,13),F(25,169)),(F(0),F(25,169)),(F(-2,17),F(64,289))]
 for batch,ns in enumerate(batches):
  for n in ns:
   for projection,perp2 in sources:
    if n==1 and perp2:continue
    assert projection*projection+perp2<F(1,4)
    square=1-perp2
    root=F(1) if not perp2 else F(12,13) if perp2==F(25,169) else F(15,17)
    assert root*root==square;t0=projection+root;start=min(F(1),t0);mass=1-start**n
    for rho in (F(0),F(1,4),F(1,2),F(3,4),F(1)):
     assert ((rho-projection)**2+perp2>=1)==(rho>=t0);checks+=1
    for threshold in (0,1,2,4,8,16,32,64):
     cut=F(1,2**threshold);tail=max(F(0),cut-start**n)
     assert tail<=cut*mass;checks+=1
     # Increasing radial masks only raise the radial lower cutoff.
     for mask_start in (F(0),F(1,4),F(3,4),F(63,64)):
      b=max(start,mask_start);masked_mass=1-b**n;masked_tail=max(F(0),cut-b**n)
      assert masked_tail<=cut*masked_mass;checks+=1
    repeat=sum((max(F(0),F(1,2**i)-start**n) for i in range(n)),F(0))
    assert repeat<=2*mass;checks+=1
    rows.append(dict(batch=batch,n=n,z_projection=str(projection),z_perpendicular_square=str(perp2),ray_outer_start=str(t0),radial_outer_mass=str(mass),canonical_repetition_conditional=str(repeat/mass) if mass else None))
 return rows,checks

def arbitrary_mask_obstruction():
 rows=[]
 for batch,ns in enumerate(((1,2,4),(8,16,32),(64,128,256))):
  for n in ns:
   rho=F(1,2)+F(1,16*n);epsilon=F(1,64*n);delta=F(1,256*n)
   xupper=rho+delta;zupper=F(1,2)-epsilon+delta;distance_lower=rho+F(1,2)-epsilon-2*delta
   assert xupper**n<F(1,2**(n-1)) and rho-delta>F(1,2)
   assert zupper<F(1,2) and distance_lower>1
   # Positive-volume x/z balls with these centers have every repetition m=n.
   rows.append(dict(batch=batch,n=n,R='1',r='1/2',gap=n,x_center_radius=str(rho),z_center_radius=str(F(1,2)-epsilon),mask_ball_radius=str(delta),x_radius_upper=str(xupper),z_radius_upper=str(zupper),outer_distance_lower=str(distance_lower),exact_masked_repetition=n,actual_E_record_realization_proved=False))
 # All-n rational scalar margin: (1+17/(128n))^n<=exp(17/128)<=128/111<2.
 assert F(128,111)<2
 return rows

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();began=time.monotonic();out=dict(status='running',scope='whole-ball geometry only; arbitrary masks kept separate from actual E records',whole_ball=[])
 with localcontext() as ctx:
  ctx.prec=120;nodes64=gauss(64);nodes96=gauss(96)
  for batch,ns in enumerate(((1,3,5),(9,17,33),(65,97,129))):
   for n in ns:
    low=whole_ball(n,64,nodes64);high=whole_ball(n,96,nodes96)
    error=max(abs(D(low[k])-D(high[k])) for k in ('outer_probability','conditional_depth_mean','conditional_repetition_partial'))
    error=max(error,max(abs(D(low['conditional_depth_tail'][t])-D(high['conditional_depth_tail'][t])) for t in low['conditional_depth_tail']))
    assert error<D('1e-15');high.update(batch=batch,cross_quadrature_absolute_error=str(error),arithmetic_precision_digits=120,interval_certificate=False);out['whole_ball'].append(high);a.output.write_text(json.dumps(out)+'\n');print(batch,n,'outside',float(D(high['outer_probability'])),'EL',float(D(high['conditional_depth_mean'])),'repeat',float(D(high['conditional_repetition_partial'])),'error',str(error),flush=True)
 out['exact_rays'],out['exact_ray_checks']=exact_rays();out['mask_obstruction']=arbitrary_mask_obstruction();out.update(status='passed',batch_counts=[3,3,3],elapsed_seconds=time.monotonic()-began,sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),whole_ball_analytic_tail='2^(-a)',whole_ball_analytic_repetition_bound=2,arbitrary_mask_uniform_repetition_bound=False,actual_E_uniform_bound_proved=False)
 a.output.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
