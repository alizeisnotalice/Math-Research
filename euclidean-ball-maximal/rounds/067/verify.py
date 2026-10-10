#!/usr/bin/env python3
"""Exact tests of exterior inversion, maximal envelopes, and coverage obstruction."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import random, json, hashlib, sys

def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def norm(a): return dot(a,a)
def determinant(a):
 a=[list(row) for row in a]; value=F(1)
 for k in range(len(a)):
  i=next((i for i in range(k,len(a)) if a[i][k]),None)
  if i is None:return F(0)
  if i!=k:a[i],a[k]=a[k],a[i];value=-value
  pivot=a[k][k];value*=pivot
  for i in range(k+1,len(a)):
   v=a[i][k]/pivot
   for j in range(k+1,len(a)):a[i][j]-=v*a[k][j]
 return value

def angular(n,a):
 # Independent direct spherical-coordinate integral, n odd.
 if not a:return F(1)
 if n==1:return ((1-a)**(-2)+(1+a)**(-2))/2
 m=(n-3)//2;p=[F(1)];base=[4*a*a-(1+a*a)**2,2*(1+a*a),F(-1)]
 for _ in range(m):
  nxt=[F(0)]*(len(p)+2)
  for i,v in enumerate(p):
   for j,w in enumerate(base):nxt[i+j]+=v*w
  p=nxt
 lo=(1-a)**2;hi=(1+a)**2
 value=sum((v*(hi**(j-n+1)-lo**(j-n+1))/F(j-n+1) for j,v in enumerate(p)),F(0))/(2*a)**(2*m+1)
 normal=sum((F(2*comb(m,j)*(-1)**j,2*j+1) for j in range(m+1)),F(0))
 return value/normal

def root_interval(value,n,bits=80):
 lo=0;hi=2**bits
 while hi-lo>1:
  mid=(lo+hi)//2
  if F(mid,2**bits)**n<=value:lo=mid
  else:hi=mid
 return F(lo,2**bits),F(hi,2**bits)

def full_line_level(points,weights,alpha):
 data=sorted(zip(points,weights));intervals=[]
 for i in range(len(data)):
  mass=F(0)
  for j in range(i,len(data)):
   mass+=data[j][1];r=mass/(2*alpha)
   lo=data[j][0]-r;hi=data[i][0]+r
   if lo<hi:intervals.append((lo,hi))
 merged=[]
 for lo,hi in sorted(intervals):
  if not merged or lo>merged[-1][1]:merged.append([lo,hi])
  else:merged[-1][1]=max(hi,merged[-1][1])
 return sum((hi-lo for lo,hi in merged),F(0))

def main():
 out=[]
 for batch,ns in enumerate(((1,2,3,4),(8,16,32),(64,128,256))):
  count=0;dets=0;maximal=0;inversions=0
  for n in ns:
   rng=random.Random(670000+n)
   for _ in range(8):
    b=tuple(F(rng.randint(-4,4),8*n) for _ in range(n));x=tuple(F(rng.randint(-16,16),4) for _ in range(n))
    v=sub(x,b);v2=norm(v)
    if not v2:continue
    z=tuple(t/v2 for t in v);A=1-norm(b);center=tuple(t/A for t in b)
    lhs=1/A**2-norm(sub(z,center));rhs=(norm(x)-1)/(A*v2)
    assert lhs==rhs
    assert tuple(t/norm(z) for t in z)==v
    t=tuple(F(rng.randint(-8,8),8) for _ in range(n));jt=tuple((t[j]-2*v[j]*dot(v,t)/v2)/v2 for j in range(n))
    assert norm(jt)==norm(t)/v2**2
    if n<=4:
     jac=[[F(int(i==j))/v2-2*v[i]*v[j]/v2**2 for j in range(n)] for i in range(n)]
     assert abs(determinant(jac))==v2**(-n);dets+=1
    inversions+=1
   # Exact finite tests of the all-radius envelope; even n avoids radicals.
   if n%2==0:
    points=[tuple(F(rng.randint(-4,4),4*n) for _ in range(n)) for _ in range(8)]
    raw=[rng.randint(1,32) for _ in points];weights=[F(w,sum(raw)) for w in raw]
    assert max(map(norm,points))<=F(1,n)
    for scale in (F(3,4),F(1),F(5,4),F(2)):
     x=(scale,)+(F(0),)*(n-1);group={}
     for y,w in zip(points,weights):
      d=norm(sub(x,y));assert d>0;group[d]=group.get(d,F(0))+w
     mass=F(0);maximum=F(0)
     for d,w in sorted(group.items()):mass+=w;maximum=max(maximum,mass/d**(n//2))
     h=sum((w/d**(n//2) for d,w in group.items()),F(0));square=sum((w/d**n for d,w in group.items()),F(0))
     assert maximum<=h and h*h<=square
     maximal+=1
   assert 4*(1+F(1,n))**n<11
   count+=1
  angular_count=0
  for n in ((1,3,5),(9,17,33),(65,97,129))[batch]:
   for a in (F(0),F(1,5),F(1,2),F(9,10)):
    assert angular(n,a)==(1+a*a)/(1-a*a)**(n+1);angular_count+=1
  line=[];N=(4,16,64)[batch]
  for trial in range(5):
   rng=random.Random(671000+100*batch+trial)
   points=[F(rng.randint(-128,128),128) for _ in range(N)]
   raw=[rng.randint(1,64) for _ in points];weights=[F(w,sum(raw)) for w in raw]
   for alpha in (F(1,8),F(1,2),F(2)):
    volume=full_line_level(points,weights,alpha);rho=max(map(abs,points));r=1/(2*alpha)
    assert (alpha*volume)**2<=4*(1+rho*rho/(r*r))
   line.append(dict(seed=671000+100*batch+trial,N=N))
  smooth=[]
  for n in ((3,4),(16,32),(128,256))[batch]:
   al,ah=root_interval(F(1,24),n);bl,bh=root_interval(F(1,8),n)
   eta=min(al,1-ah)/8;assert eta>0 and ah+2*eta<1 and bh+2*eta<2
   D=n+6
   assert F(2)**(1-D)<F(1,24)/2**n
   assert 2*eta<al/2 and n**n>F(9,4)
   smooth.append(dict(n=n,D=D,eta_over_R1=str(eta),g1_over_alpha='1/16',plateau_over_alpha='3/2'))
  out.append(dict(batch=batch,dimensions=list(ns),inversion_checks=inversions,determinants=dets,all_radius_envelope_checks=maximal,angular_integrals=angular_count,line_inputs=line,line_threshold_checks=15,smooth_obstructions=smooth))
 path=Path(__file__)
 print(json.dumps(dict(round=67,status='passed',python=sys.version.split()[0],batches=out,sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest()},main_theorem_proved=False,scope='Exact finite implementation checks; general estimates and obstructions are proved in report.md.'),indent=2))
if __name__=='__main__':main()
