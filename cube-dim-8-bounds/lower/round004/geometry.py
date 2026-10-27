#!/usr/bin/env python3
from pathlib import Path
from itertools import product
from math import prod, factorial
from fractions import Fraction as F
import hashlib,json,time,resource,sys
H=Path(__file__).resolve().parent
D=10**35; S=4*D; LIMIT=2**20
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def root(n):
 a,b=0,1
 while b**8<=n:b*=2
 while b-a>1:
  c=(a+b)//2
  if c**8<=n:a=c
  else:b=c
 assert a**8<=n<(a+1)**8
 return a
class Exhausted(Exception):pass
class Union:
 def __init__(self,iv,budget):self.iv=iv;self.cache={};self.budget=budget
 def call(self,state):
  if not state:return 0
  if state in self.cache:return self.cache[state]
  if self.budget[0]>=LIMIT:raise Exhausted
  self.budget[0]+=1
  if len(state[0])==1:
   ints=sorted(self.iv[b[0]] for b in state);a,b=ints[0];total=0
   for c,d in ints[1:]:
    if b<c:total+=b-a;a,b=c,d
    else:b=max(b,d)
   total+=b-a
  else:
   labels={t for box in state for t in box};cuts=sorted({z for t in labels for z in self.iv[t]});total=0
   for a,b in zip(cuts,cuts[1:]):
    child=set()
    for box in state:
     for t in set(box):
      l,h=self.iv[t]
      if l<=a and b<=h:
       z=list(box);z.remove(t);child.add(tuple(z))
    if child:total+=(b-a)*self.call(tuple(sorted(child)))
  self.cache[state]=total;return total

def build():
 masses={prod(e) for e in product([1,2],*[range(1,5)]*7)}
 rr={m:root(m*D**8) for m in masses};tail=(2*(D+rr[2]),2*(D+rr[2]+1));records=[]
 for e in product([1,2],*[range(1,5)]*7):
  m=prod(e);rad=(2*rr[m],2*(rr[m]+1));inner=[];outer=[]
  for ell in e[1:]:
   q=-(ell//2) if ell%2 else 1-ell//2;c=q+ell-1
   l=(max(0,c*tail[0]-rad[1]),max(0,c*tail[1]-rad[0]))
   h=(min(tail[0]//2,q*tail[1]+rad[0]),min(tail[1]//2,q*tail[0]+rad[1]))
   inner.append((l[1],h[0]));outer.append((l[0],h[1]))
  records.append(dict(ell=e,rad=rad,inner=tuple(inner),outer=tuple(outer)))
 orbit=[r for r in records if tuple(sorted(r['ell'][1:]))==r['ell'][1:]]
 assert len(records)==32768 and len(orbit)==240
 assert sum(factorial(7)//prod(factorial(r['ell'][1:].count(v)) for v in set(r['ell'][1:])) for r in orbit)==32768
 return rr,tail,records,orbit

def first(ell,aa,rad,mode):
 low=(0,0) if ell==1 else (max(0,aa[0]-rad[1]),max(0,aa[1]-rad[0]))
 high=(min(aa[0]//2,rad[0]),min(aa[1]//2,rad[1]))
 return (low[1],high[0]) if mode=='inner' else (low[0],high[1])
def calculate(event,spacing,rr,tail,records,orbit,ref,inherited):
 a=spacing
 cell=F(a[1]*tail[1]**7,S**8)
 r={'event':event,'spacing_endpoints':list(map(str,spacing)),'cell_upper':str(cell),'cell_prefilter_excluded':cell<=inherited}
 budget=[0];ans={}; t=time.time()
 for mode in ('inner','outer'):
  iv=sorted({z for b in orbit for z in b[mode]});ids={z:i for i,z in enumerate(iv)};eng=Union(iv,budget);boxes=[]
  for b in orbit:
   f=first(b['ell'][0],a,b['rad'],mode)
   if f[0]<f[1] and all(l<h for l,h in b[mode]):boxes.append((f,tuple(sorted(ids[z] for z in b[mode]))))
  cuts=sorted({z for pair,_ in boxes for z in pair});total=0
  for lo,hi in zip(cuts,cuts[1:]):
   state=tuple(sorted({box for (l,h),box in boxes if l<=lo and hi<=h}))
   total+=(hi-lo)*eng.call(state)
  ans[mode]=F(256*total,S**8)
  r[mode+'_states']=len(eng.cache)
 assert ans['inner']<=ans['outer']<=cell
 r.update(union_interval=[str(ans['inner']),str(ans['outer'])],states=budget[0],wall_seconds=time.time()-t,peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 if True:
  positive=[];checks=0;width=[];outer=0
  for b in records:
   ins=(first(b['ell'][0],a,b['rad'],'inner'),)+b['inner'];outs=(first(b['ell'][0],a,b['rad'],'outer'),)+b['outer']
   if all(l<h for l,h in outs):outer+=1
   if not all(l<h for l,h in ins):continue
   positive.append(b)
   for i,(ell,(l,h)) in enumerate(zip(b['ell'],ins)):
    spacing=a if i==0 else tail;q=-(ell//2) if ell%2 else 1-ell//2
    assert l>=(q+ell-1)*spacing[1]-b['rad'][0] and h<=q*spacing[1]+b['rad'][0]
    width.append(h-l);checks+=2
  eps=F(1,10**25);theta=1-F(1,10**18);K=10**19
  assert F(min(width),S)>2*eps and max(a[1],tail[1])<4*S
  loss=256*len(positive)*16*eps*2**7
  finite=theta*(ans['inner']-loss)*F(2*K+1,2*K+5)**8
  r.update(positive_inward_raw_boxes=len(positive),positive_outward_raw_boxes=outer,containment_checks=checks,minimum_width=str(F(min(width),S)),finite_L1=dict(epsilon=str(eps),eta=str(eps/2),theta=str(theta),K=K,collar=2,loss_upper=str(loss),lower=str(finite),margin=str(finite-ref)))
  r['winner']=finite>ref
 return r

