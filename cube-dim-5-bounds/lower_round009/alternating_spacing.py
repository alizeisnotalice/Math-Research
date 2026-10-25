#!/usr/bin/env python3
"""Exact finite period-two first-axis spacing certificate, with unit weights."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
LIMIT_SECONDS,LIMIT_NODES=120,200000
def lift(F,q):return F(q.numerator)/q.denominator
def root(n,u,v):
 e2=e3=0
 while n%2==0:n//=2;e2+=1
 while n%3==0:n//=3;e3+=1
 assert n==1;return u**e2*v**e3
def decode(F,d):
 a=[Q(0)]*25
 for k,x in d.items():i,j=map(int,k.split(','));a[5*i+j]=Q(x)
 return F(a=a)
def run(delta):
 started=time.monotonic();F,u,v,bounds,_,_=field();a1=1+u**3/2;s=(1+u)/2;deltaF=lift(F,delta)
 def atom(j): return j*a1+(deltaF if j%2 else F(0))
 raw=[]
 for ls in product(range(1,5),repeat=5):
  r=root(prod(ls),u,v)/2;l1=ls[0]
  for j in range(-6,6):
   lo=max(F(0),atom(j+l1-1)-r);hi=min(2*a1,atom(j)+r)
   if hi<=lo:continue
   box=[lo,hi]
   for l in ls[1:]:
    if l%2:x,y=F(0),min(s/2,r-(l-1)*s/2)
    else:x,y=max(F(0),l*s/2-r),s/2
    if y<=x:break
    box.extend((x,y))
   else:raw.append((j,ls,tuple(box)))
 unique={b:(j,ls) for j,ls,b in raw};kept=[]
 for b,val in unique.items():
  if not any(o!=b and all(o[2*k]<=b[2*k] and b[2*k+1]<=o[2*k+1] for k in range(5)) for o in unique):kept.append((val,b))
 retained=tuple(sorted(b for _,b in kept));nodes=0
 @lru_cache(None)
 def rec(boxes):
  nonlocal nodes
  nodes+=1
  if nodes>LIMIT_NODES:raise RuntimeError('node cap')
  if time.monotonic()-started>LIMIT_SECONDS:raise TimeoutError('wall cap')
  if not boxes:return F(0)
  if len(boxes[0])==2:
   out=F(0);lo=hi=None
   for x,y in sorted(boxes):
    if lo is None:lo,hi=x,y
    elif hi<x:out+=hi-lo;lo,hi=x,y
    elif hi<y:hi=y
   return out+hi-lo
  cuts=sorted({x for b in boxes for x in b[:2]});out=F(0)
  for lo,hi in zip(cuts,cuts[1:]):out+=(hi-lo)*rec(tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]})))
  return out
 volume=rec(retained);segments=[]
 for k in range(5):
  cuts=sorted({b[2*k+q] for b in retained for q in (0,1)});one=[]
  for lo,hi in zip(cuts,cuts[1:]):
   mask=sum(1<<i for i,b in enumerate(retained) if b[2*k]<=lo and hi<=b[2*k+1])
   if mask:one.append((hi-lo,mask))
  segments.append(one)
 @lru_cache(None)
 def grid(k,mask):
  if k==5:return F(1)
  out=F(0)
  for d,nxt in segments[k]:
   both=nxt if mask<0 else mask&nxt
   if both:out+=d*grid(k+1,both)
  return out
 assert volume==grid(0,-1); ratio=8*volume
 c=json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text());C=decode(F,c['coefficients']);diff=ratio-C
 lo,hi=bounds(ratio.a);dl,dh=bounds(diff.a)
 return {'status':'PASS_exact_improvement' if dl>0 else 'NO_IMPROVEMENT_exact','delta':str(delta),'first_axis_period_atoms':['0','a1+delta'],'first_axis_start_range':[-6,5],
  'raw_boxes':len(raw),'retained_boxes':len(kept),'strict_level':'r^5=prod(lengths)/32 (open boundary handled by threshold limit)',
  'receiver_domain':'[0,2*a1] x prod tail [0,s/2]','period_mass':2,'ratio_normalization':'8*absolute_union_volume','union_engines_agree':True,
  'coordinate_segment_counts':[len(x) for x in segments],'recursive_nodes':nodes,'wall_seconds':time.monotonic()-started,
  'ratio_coefficients':{f'{i//5},{i%5}':str(x) for i,x in enumerate(ratio.a) if x},'ratio_interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)],
  'difference_from_exact_C_interval':[str(dl),str(dh)],'strictly_above_exact_C':dl>0,'C_source_sha256':hashlib.sha256((ROOT/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest(),
  'limits':{'wall_seconds':LIMIT_SECONDS,'recursive_nodes':LIMIT_NODES},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
def main():
 if len(sys.argv)!=2:raise SystemExit('usage: alternating_spacing.py delta')
 d=Q(sys.argv[1]);assert d in (Q(1,32),Q(1,8),Q(1,4))
 try:out=run(d)
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','delta':str(d),'reason':str(e)}
 (HERE/f'alternating_spacing_delta_{d.numerator}_{d.denominator}_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
