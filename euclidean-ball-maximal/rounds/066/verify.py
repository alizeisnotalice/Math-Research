#!/usr/bin/env python3
"""Rational enclosing-ball certificates and geometry-weighted shell budgets."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb
import importlib.util,json,hashlib,random,sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
def module(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def norm(a):return dot(a,a)
def solve(a,b):
 n=len(b);a=[list(row)+[v] for row,v in zip(a,b)]
 for k in range(n):
  pivot=next((i for i in range(k,n) if a[i][k]),None)
  if pivot is None:return None
  a[k],a[pivot]=a[pivot],a[k];v=a[k][k];a[k]=[x/v for x in a[k]]
  for i in range(n):
   if i!=k:
    v=a[i][k];a[i]=[x-v*y for x,y in zip(a[i],a[k])]
 return [a[i][-1] for i in range(n)]
def minimum_ball(points):
 n=len(points[0])
 for k in range(1,min(n+1,len(points))+1):
  for ids in combinations(range(len(points)),k):
   a=points[ids[0]];vs=[sub(points[i],a) for i in ids[1:]]
   beta=solve([[dot(v,z) for z in vs] for v in vs],[norm(v)/2 for v in vs])
   if beta is None:continue
   weights=[1-sum(beta)]+beta
   if min(weights)<0:continue
   c=tuple(a[j]+sum((v[j]*b for v,b in zip(vs,beta)),F(0)) for j in range(n));rho=norm(sub(c,a))
   if all(norm(sub(x,c))<=rho for x in points):return c,rho,ids,weights
 raise AssertionError('no enclosing-ball certificate')
def check_certificate(points,c,rho,ids,weights,seed):
 n=len(c);assert len(ids)<=n+1 and sum(weights)==1 and min(weights)>=0
 assert all(sum((w*points[i][j] for i,w in zip(ids,weights)),F(0))==c[j] for j in range(n))
 assert all(norm(sub(points[i],c))==rho for i in ids)
 assert all(norm(sub(a,c))<=rho for a in points)
 rng=random.Random(seed);support=points if len(points)<=10 else points[:8]+[points[-1]]
 raw=[rng.randint(1,32) for _ in support];prob=[F(v,sum(raw)) for v in raw]
 mean=tuple(sum((w*a[j] for w,a in zip(prob,support)),F(0)) for j in range(n))
 variance=sum((w*norm(sub(a,mean)) for w,a in zip(prob,support)),F(0))
 pair=sum((v*w*norm(sub(a,b))/2 for v,a in zip(prob,support) for w,b in zip(prob,support)),F(0))
 assert pair==variance and variance+norm(sub(mean,c))<=rho
 for _ in range(8):
  x=tuple(F(rng.randint(-16,16),4) for _ in range(n));far=max(norm(sub(x,a)) for a in points)
  assert norm(sub(x,c))+rho<=far
  assert norm(sub(x,mean))+variance<=far
 return dict(n=n,points=len(points),contact_count=len(ids),radius_squared=str(rho),variance=str(variance),seed=seed)
def shell_bound(n,r,R,a):
 lo=max(r,a)
 if lo>=R:return F(0)
 m=(n-1)//2
 def prim(t):return sum((comb(m,j)*(-a*a)**(m-j)*t**(2*j+1)/F(2*j+1) for j in range(m+1)),F(0))
 return n*(prim(R)-prim(lo))
def outward(value,up=False,bits=80):
 scale=2**bits;num=value.numerator*scale;den=value.denominator
 return F((num+den-1)//den if up else num//den,scale)
def main():
 old=module('rounds/064/verify.py','old66');certs=[];shells=[];preload=[]
 for n in (1,2,3,4):
  for serial in range(3):
   seed=660000+100*n+serial;rng=random.Random(seed)
   points=list(dict.fromkeys(tuple(F(rng.randint(-8,8),4) for _ in range(n)) for _ in range(n+4)))
   c,rho,ids,weights=minimum_ball(points);certs.append(dict(batch=0,**check_certificate(points,c,rho,ids,weights,seed)))
 for batch,ns in enumerate(((8,16,32),(64,128,256)),1):
  for n in ns:
   points=[tuple(F(int(i==j)) for j in range(n)) for i in range(n)]+[tuple(F(-1) for _ in range(n))]
   t=F(1-n,2*(n+1));c=(t,)*n;b=F(n+3,2*(n+1)**2);weights=[b]*n+[1-n*b]
   rho=norm(sub(points[0],c));certs.append(dict(batch=batch,**check_certificate(points,c,rho,tuple(range(n+1)),weights,660000+n)))
 for batch,ns in enumerate(((1,3,5),(9,17,33),(65,97,129))):
  for n in ns:
   checks=0;best=F(0);strongest=F(1)
   for a in (F(0),F(1,4),F(3,4),F(99,100),F(11,10)):
    for r in (F(1,4),F(1,2),F(9,10),F(999,1000)):
     actual=old.lens(n,F(1),a)-old.lens(n,r,a)
     assert old.lens(n,F(1),a)**2<=max(F(0),1-a*a)**n
     for sigma in (a,a/2):
      budget=shell_bound(n,r,F(1),sigma)
      theta=(1-sigma*sigma)**((n-1)//2) if sigma<1 else F(0)
      assert 0<=actual<=budget<=(1-r**n)*theta
      if budget:best=max(best,actual/budget)
      if a<1 and sigma==a:strongest=min(strongest,theta)
      checks+=1
   shells.append(dict(batch=batch,n=n,checks=checks,max_actual_over_integral=str(best),smallest_active_theta=str(strongest)))
 for batch,ms in enumerate(((3,4),(5,6),(7,8))):
  for m in ms:
   n=2**m;eta=F(1,2**(20+3*m));cl,ch=old.root_interval(F(3,2),n);rl,rh=old.root_interval(F(2),n)
   gamma_lo=cl+eta;gamma_hi=ch+eta
   assert gamma_hi<rl and gamma_lo**n>F(3,2)
   preload_lo=gamma_lo**n/2;preload_hi=gamma_hi**n/2
   assert F(3,4)<preload_lo<=preload_hi<1
   factor_squared_upper=1-(gamma_lo/rh)**2
   assert 0<factor_squared_upper<F(1,n)
   # theta is a half-integral power; certify the squared final coefficient.
   squared_cap=F(3,16)**2*(n-1)**2*factor_squared_upper**(n-1)
   analytic_squared=F(3,16)**2*F((n-1)**2,n**(n-1))
   assert squared_cap<analytic_squared<1
   preload.append(dict(batch=batch,n=n,eta_over_R=str(eta),preload_fraction_lower=str(outward(preload_lo)),preload_fraction_upper=str(outward(preload_hi,True)),theta_squared_upper_base=str(factor_squared_upper),capacity_bound_squared=str(analytic_squared)))
 deps=json.loads((ROOT/'rounds/065/verification.json').read_text())['sha256'];deps['rounds/066/verify.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 for p,h in deps.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 print(json.dumps(dict(round=66,status='passed',python=sys.version.split()[0],contact_certificates=certs,shells=shells,preload=preload,sha256=deps,main_theorem_proved=False,scope='Exact finite geometric certificates; no global weighted-tree covering theorem is asserted.'),indent=2))
if __name__=='__main__':main()
