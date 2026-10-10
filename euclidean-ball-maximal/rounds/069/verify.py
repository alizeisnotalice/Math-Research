#!/usr/bin/env python3
"""Exact truncation algebra and complete one-dimensional step-density levels."""
from fractions import Fraction as F
from pathlib import Path
import random,json,hashlib,sys

def mass(data):return sum(((b-a)*h for a,b,h in data),F(0))
def primitive(data,x):return sum((h*max(F(0),min(x,b)-a) for a,b,h in data),F(0))
def merge(items):
 out=[]
 for a,b in sorted(items):
  if a>=b:continue
  if out and a<=out[-1][1]:out[-1][1]=max(b,out[-1][1])
  else:out.append([a,b])
 return out

def length(items):return sum((b-a for a,b in merge(items)),F(0))
def levels(data,alpha):
 data=[v for v in data if v[2]>0]
 if not data:return []
 ends=sorted({x for a,b,h in data for x in (a,b)})
 reach=mass(data)/(2*alpha)
 breaks=sorted({ends[0]-reach,ends[-1]+reach,*ends,*((a+b)/2 for a in ends for b in ends)})
 out=[]
 for l,r in zip(breaks,breaks[1:]):
  mid=(l+r)/2
  if any(a<mid<b and h>alpha for a,b,h in data):out.append((l,r));continue
  for e in ends:
   def excess(x):
    d=abs(x-e)
    return primitive(data,x+d)-primitive(data,x-d)-2*alpha*d
   a=excess(l);b=excess(r)
   if a>=0 and b>=0 and a+b>0:out.append((l,r))
   elif a>0 or b>0:
    root=l-a*(r-l)/(b-a)
    out.append((l,root) if a>0 else (root,r))
 return merge(out)

def main():
 batches=[]
 for batch,powers in enumerate(((1,2,3),(4,6,8),(12,20,32))):
  algebra=0;counts=0;changed=0
  for k in powers:
   delta=F(1,2**k);a=delta*delta/3;tau=1-delta
   for eta in (F(0),F(1,1000),F(1,10),F(2,3)):
    bound=3*(delta+eta*tau)/(2+delta)
    assert bound<=F(3,2)*(delta+eta)
    for u in (F(0),min(bound,F(9,10))/2,min(bound,F(9,10))):
     assert (1-eta)<=a*u/delta**2+(1-u)/tau
     assert u*(1-a*tau/delta**2)<=delta+eta*tau
     residual=tau*(1-eta-u/3)/(1-u)
     assert residual>=tau*(1-eta)
     algebra+=1
  N=(1,2,4)[batch]
  for trial in range(3):
   seed=690000+100*batch+trial;rng=random.Random(seed)
   data=[]
   for j in range(N):
    left=F(3*j);width=F(rng.randint(1,4),8);height=F(rng.randint(1,16))
    data.append((left,left+width,height))
   W=mass(data);data=[(l,r,h/W) for l,r,h in data]
   for alpha in (F(1,8),F(1,2),F(2)):
    original=levels(data,alpha);v=length(original)
    for L in (F(1,2),F(2)):
     g=[(l,r,min(h,L*alpha)) for l,r,h in data]
     hdata=[(l,r,max(F(0),h-L*alpha)) for l,r,h in data]
     u=mass(g);assert u+mass(hdata)==1
     carrier=sum((r-l for l,r,h in data if h>L*alpha),F(0))
     assert L*alpha*carrier<=u
     for delta in (F(1,3),F(2,3)):
      low=levels(g,delta*alpha);high=levels(hdata,(1-delta)*alpha)
      cover=merge(low+high)
      assert length(cover+original)==length(cover)
      assert v<=length(low)+length(high)
      if u<1:
       normalized=[(l,r,h/(1-u)) for l,r,h in hdata]
       assert levels(normalized,(1-delta)*alpha/(1-u))==high
       changed+=int(high!=original)
      counts+=1
  # Analytic single-interval benchmark, independent of the cell subdivision.
  for alpha in (F(1,4),F(1,2),F(3,4),F(1),F(2)):
   got=levels([(F(0),F(1),F(1))],alpha)
   expected=[] if alpha>=1 else [[1-1/(2*alpha),1/(2*alpha)]] if alpha<=F(1,2) else [[F(0),F(1)]]
   # For alpha>1/2 the exterior halo disappears, while the interior remains.
   assert got==expected
  batches.append(dict(batch=batch,delta_powers=list(powers),algebra_checks=algebra,step_intervals=N,seeds=[690000+100*batch+t for t in range(3)],complete_level_split_checks=counts,changed_residual_levels=changed,single_interval_checks=5))
 path=Path(__file__)
 print(json.dumps(dict(round=69,status='passed',python=sys.version.split()[0],batches=batches,sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest()},main_theorem_proved=False,scope='Conditional near-extremizer algebra and exact one-dimensional truncation; no divergent near-extremizer constructed.'),indent=2))
if __name__=='__main__':main()
