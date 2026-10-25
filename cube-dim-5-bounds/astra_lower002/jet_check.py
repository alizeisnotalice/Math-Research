#!/usr/bin/env python3
"""Independent mixed-capture jet union at the registered cell midpoint.

Imports the established scalar Q(u,v) field, not exposed-facet code.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import json,hashlib,time,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
F,u,v,bounds,_,_=field();ZERO=F(0)
def decode(o):
 a=[Q(0)]*25
 for key,value in o.items():i,j=map(int,key.split(','));a[5*i+j]=Q(value)
 return F(a=a)
def encode(x):return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.a) if c}
class Jet:
 __slots__=('value','speed')
 def __init__(self,value=0,speed=0):self.value=value if isinstance(value,F) else F(value);self.speed=speed if isinstance(speed,F) else F(speed)
 def __hash__(self):return hash((self.value,self.speed))
 def __eq__(self,o):return self.value==o.value and self.speed==o.speed
 def __lt__(self,o):return self.value<o.value if self.value!=o.value else self.speed<o.speed
 def __le__(self,o):return self==o or self<o
 def __add__(self,o):return Jet(self.value+o.value,self.speed+o.speed)
 def __sub__(self,o):return Jet(self.value-o.value,self.speed-o.speed)
 def __mul__(self,o):return Jet(self.value*o.value,self.value*o.speed+self.speed*o.value)
JZERO=Jet();JONE=Jet(1)
A=(1+u**3/2,)+((1+u)/2,)*4
reg=json.loads((HERE/'event_cell.json').read_text());T=decode(reg['nodes'][1]['position_coefficients'])
@lru_cache(None)
def radius(n,m):
 t=n;i=j=0
 while t%2==0:t//=2;i+=1
 while t%3==0:t//=3;j+=1
 assert t==1
 r=u**i*v**j/2
 return Jet(r,r*m/(5*n))
def closest(axis,shift,x,n):
 # Select actual coordinate indices by exact distance from a segment interior.
 def distance(k):
  d=k*A[axis]+shift-x
  return -d if d<ZERO else d
 return tuple(sorted(sorted(range(-4,5),key=distance)[:n]))
begin=time.monotonic();nodes=0
def budget():
 if time.monotonic()-begin>120 or nodes>200000:raise TimeoutError('registered jet cap')
cuts=(ZERO,T,A[1]/2,T+A[1]/2,A[1]);raw={};configs=0
for lo,hi in zip(cuts,cuts[1:]):
 x=(lo+hi)/2
 sets={}
 for axis in range(5):
  for orbit in (0,1):
   shift=T if axis==1 and orbit else ZERO
   point=x if axis==1 else A[axis]/4
   for n in range(1,5):sets[axis,orbit,n]=closest(axis,shift,point,n)
 def add(l0,l1):
  global configs
  configs+=1
  n=prod(l0);m=prod(l1) if l1 else 0;r=radius(n,m);b=[]
  for axis in range(5):
   allcoords=[]
   for orb,ls in ((0,l0),(1,l1)):
    if ls:
     shift=T if orb and axis==1 else ZERO
     allcoords.extend(k*A[axis]+shift for k in sets[axis,orb,ls[axis]])
   lower=max(Jet(lo if axis==1 else ZERO),Jet(max(allcoords))-r)
   upper=min(Jet(hi if axis==1 else A[axis]/2),Jet(min(allcoords))+r)
   if upper<=lower:return
   b.extend((lower,upper))
  raw[tuple(b)]={'l0':l0,'l1':l1,'n':n,'m':m,'sets':{(axis,orb):sets[axis,orb,ls[axis]] for orb,ls in ((0,l0),(1,l1)) if ls for axis in range(5)}}
 for ls in product(range(1,5),repeat=5):add(ls,None)
 opts=[[(p,q) for p in range(1,5) for q in range(1,5) if (abs(p-q)<=1 if axis==1 else p==q)] for axis in range(5)]
 for pq in product(*opts):add(tuple(y[0] for y in pq),tuple(y[1] for y in pq))
 budget()
kept=tuple(sorted(b for b in raw if not any(o!=b and all(o[2*j]<=b[2*j] and b[2*j+1]<=o[2*j+1] for j in range(5)) for o in raw)))
guard_count=0;selected_atoms=0
assert all(T!=k*A[1] for k in range(-8,9))
for b in kept:
 data=raw[b];r=radius(data['n'],data['m']);atom_tags=[]
 for orb,key in ((0,'l0'),(1,'l1')):
  ls=data[key]
  if ls is None:continue
  tags=list(product(*(data['sets'][axis,orb] for axis in range(5))))
  assert len(tags)==prod(ls)==len(set(tags));atom_tags.extend((orb,tag) for tag in tags)
  for axis in range(5):
   shift=T if orb and axis==1 else ZERO
   for k in data['sets'][axis,orb]:
    p=Jet(k*A[axis]+shift)
    assert p-r<=b[2*axis] and b[2*axis+1]<=p+r
    guard_count+=2
 assert len(set(atom_tags))==data['n']+data['m'];selected_atoms+=len(atom_tags)
@lru_cache(None)
def volume(boxes):
 global nodes
 nodes+=1;budget()
 if not boxes:return JZERO
 cuts=sorted({x for b in boxes for x in b[:2]});total=JZERO
 for lo,hi in zip(cuts,cuts[1:]):
  sub=tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]}))
  if not sub:continue
  total=total+(hi-lo)*(JONE if not sub[0] else volume(sub))
 return total
V=volume(kept);C=decode(json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients'])
assert 16*V.value==C
D=16*V.speed-C
facet=decode(json.loads((HERE/'cell_derivative_results.json').read_text())['node_derivatives'][1]['coefficients'])
assert D==facet
out={'status':'EXACT_MIXED_BOX_JET_AGREES_WITH_EXPOSED_FACET_KERNEL','node_index':1,'normalization':'16*half-reflected-cell volume/(1+epsilon)','constant_exactly_C':True,'derivative_coefficients':encode(D),'derivative_interval':[str(x) for x in bounds(D.a)],'derivative_display_only':[float(x) for x in bounds(D.a)],'configs':configs,'raw_jet_boxes':len(raw),'retained_jet_boxes':len(kept),'selected_atom_tags_checked':selected_atoms,'coordinate_atom_endpoint_jet_inequalities':guard_count,'distinct_orbits_proved_by_nonintegral_shift':True,'nodes':nodes,'seconds':time.monotonic()-begin,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_imports':'scalar field only; no exposed-facet or union producer imports'}
(HERE/'jet_check_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('derivative_coefficients','derivative_interval')},indent=2))
