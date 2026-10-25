#!/usr/bin/env python3
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import json,time,sys,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field

def calculate(spec):
 started=time.monotonic();F,u,v,bounds,_,_=field();axes=tuple(spec['axes']);nodes=0
 def cap():
  if time.monotonic()-started>120:raise TimeoutError('wall cap')
 class J:
  __slots__=('c','d','g')
  def __init__(self,c=0,d=0,g=None):self.c=c if isinstance(c,F) else F(c);self.d=d if isinstance(d,F) else F(d);self.g=tuple(F(0) for _ in range(5)) if g is None else g
  def __hash__(self):return hash((self.c,self.d,self.g))
  def __eq__(self,o):o=o if isinstance(o,J) else J(o);return self.c==o.c and self.d==o.d and self.g==o.g
  def __add__(self,o):o=o if isinstance(o,J) else J(o);return J(self.c+o.c,self.d+o.d,tuple(x+y for x,y in zip(self.g,o.g)))
  __radd__=__add__
  def __neg__(self):return J(-self.c,-self.d,tuple(-x for x in self.g))
  def __sub__(self,o):return self+-o
  def __rsub__(self,o):return J(o)+-self
  def __mul__(self,o):o=o if isinstance(o,J) else J(o);return J(self.c*o.c,self.c*o.d+self.d*o.c,tuple(self.c*y+x*o.c for x,y in zip(self.g,o.g)))
  __rmul__=__mul__
  def __truediv__(self,o):return J(self.c/o,self.d/o,tuple(x/o for x in self.g))
  def __lt__(self,o):o=o if isinstance(o,J) else J(o);return self.c<o.c if self.c!=o.c else self.d<o.d
  def __le__(self,o):return self==o or self<o
 a=(J(1+u**3/2),)+(J((1+u)/2),)*4
 @lru_cache(None)
 def radius(n,alt):
  t=n;e2=e3=0
  while t%2==0:t//=2;e2+=1
  while t%3==0:t//=3;e3+=1
  assert t==1
  r=u**e2*v**e3/2
  return J(r,r*alt/(5*n))
 @lru_cache(None)
 def interval(k,l,j,n,alt):
  r=radius(n,alt);aa=a[k]
  lo,hi=max(J(0),(j+l-1)*aa-r),min(aa if k in axes else aa/2,j*aa+r)
  def mark(x):return J(x.c,x.d,tuple(x.d if i==k else F(0) for i in range(5)))
  return mark(lo),mark(hi)
 raw=[];masschecks=0;zero_constant_boxes=0
 for ls in product(range(1,5),repeat=5):
  cap();n=prod(ls)
  for js in product(range(-3,4),repeat=2):
   alt=((-1 if sum(js)%2 else 1)*prod(ls[k] for k in range(5) if k not in axes)) if all(ls[k]%2 for k in axes) else 0
   box=[];starts=[]
   for k,l in enumerate(ls):
    j=js[axes.index(k)] if k in axes else (-(l-1)//2 if l%2 else 1-l//2)
    lo,hi=interval(k,l,j,n,alt)
    if hi<=lo:break
    box.extend((lo,hi));starts.append(j)
   else:
    direct=sum((-1 if (x+y)%2 else 1) for x,y in product(range(starts[axes[0]],starts[axes[0]]+ls[axes[0]]),range(starts[axes[1]],starts[axes[1]]+ls[axes[1]])))*prod(ls[k] for k in range(5) if k not in axes)
    assert direct==alt;masschecks+=1
    if any(box[2*k].c==box[2*k+1].c for k in range(5)):zero_constant_boxes+=1
    raw.append((tuple(box),{'lengths':ls,'starts':starts,'distinct_atom_count':n,'mass_derivative':alt}))
 unique={b:m for b,m in raw};kept=[]
 for b in unique:
  cap()
  if not any(o!=b and all(o[2*k]<=b[2*k] and b[2*k+1]<=o[2*k+1] for k in range(5)) for o in unique):kept.append(b)
 @lru_cache(None)
 def sweep(boxes):
  nonlocal nodes
  nodes+=1;cap()
  if nodes>200000:raise RuntimeError('node cap')
  if not boxes:return J(0)
  if len(boxes[0])==2:
   total=J(0);lo=hi=None
   for x,y in sorted(boxes):
    if lo is None:lo,hi=x,y
    elif x>hi:total+=hi-lo;lo,hi=x,y
    else:hi=max(hi,y)
   return total+hi-lo
  cuts=sorted({x for b in boxes for x in b[:2]});total=J(0)
  for lo,hi in zip(cuts,cuts[1:]):total+=(hi-lo)*sweep(tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]})))
  return total
 vol=sweep(tuple(sorted(kept)));segments=[];events=[]
 for k in range(5):
  cuts=sorted({b[2*k+d] for b in kept for d in (0,1)});ss=[]
  for lo,hi in zip(cuts,cuts[1:]):
   mask=sum(1<<i for i,b in enumerate(kept) if b[2*k]<=lo and hi<=b[2*k+1])
   if mask:ss.append((hi-lo,mask))
  segments.append(ss)
  groups={}
  for x in cuts:groups.setdefault(x.c,set()).add(x.d)
  events.append({'coordinate':k,'distinct_endpoint_values':len(groups),'endpoint_jet_count':len(cuts),'coincident_values_with_multiple_speeds':sum(len(z)>1 for z in groups.values())})
 @lru_cache(None)
 def masks(k,mask):
  cap()
  if k==5:return J(1)
  return sum((d*masks(k+1,mask&nxt) for d,nxt in segments[k] if mask&nxt),J(0))
 assert vol==masks(0,(1<<len(kept))-1);ratio=8*vol
 c=json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text());ca=[Q(0)]*25
 for key,val in c['coefficients'].items():i,j=map(int,key.split(','));ca[5*i+j]=Q(val)
 assert ratio.c==F(a=ca)
 dlo,dhi=bounds(ratio.d.a)
 assert sum(ratio.g,F(0))==ratio.d
 def encode(x):return {f'{i//5},{i%5}':str(v) for i,v in enumerate(x.a) if v}
 return {'status':'STRICT_NEGATIVE_RIGHT_DERIVATIVE' if dhi<0 else 'POSITIVE_RIGHT_DERIVATIVE' if dlo>0 else 'ZERO_OR_UNRESOLVED','specification':spec,'period_mass':'4','period_mass_derivative':'0','normalization':'8*volume; no derivative of denominator','constant_exactly_equals_full_C':True,'right_derivative_coefficients':encode(ratio.d),'derivative_interval':[str(dlo),str(dhi)],'derivative_display_only':[float(dlo),float(dhi)],'coordinate_face_contributions':[{'coordinate':k,'coefficients':encode(x),'interval':list(map(str,bounds(x.a))),'display_only':[float(t) for t in bounds(x.a)]} for k,x in enumerate(ratio.g)],'raw_boxes':len(raw),'retained_boxes':len(kept),'new_zero_constant_boxes':zero_constant_boxes,'direct_distinct_mass_derivative_checks':masschecks,'union_engines_agree':True,'recursive_nodes':nodes,'endpoint_events':events,'segments':[len(s) for s in segments],'seconds':time.monotonic()-started,'retained_jets':[dict(unique[b],endpoint_jets=[{'value':encode(x.c),'speed':encode(x.d)} for x in b]) for b in kept],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'C_source_sha256':hashlib.sha256((ROOT/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest()}
if __name__=='__main__':
 specs=json.loads((HERE/'registration.json').read_text())['inputs'];spec=next(s for s in specs if s['name']==sys.argv[1])
 try:out=calculate(spec)
 except (TimeoutError,RuntimeError) as e:out={'status':'COMPUTE_LIMIT','reason':str(e),'specification':spec}
 (HERE/(spec['name']+'_derivative.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['retained_jets','derivative_interval']},indent=2))
