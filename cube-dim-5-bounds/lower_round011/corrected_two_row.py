#!/usr/bin/env python3
"""Corrected two-row finite witnesses; explicit atom inclusion guards."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
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
def run(tau):
 started=time.monotonic(); F,u,v,bounds,_,_=field();a1=1+u**3/2;s=(1+u)/2;h=a1*lift(F,tau)
 raw=[]
 def budget():
  if time.monotonic()-started>LIMIT_SECONDS:raise TimeoutError("wall cap")
 for ls in product(range(1,5),repeat=5):
  budget();l1,l2=ls[:2];r=root(prod(ls),u,v)/2
  for j in range(-6,7):
   for k in range(-6,7):
    shifts=[h if q%2 else F(0) for q in range(k,k+l2)]
    minx=j*a1+min(shifts);maxx=(j+l1-1)*a1+max(shifts)
    lo=max(F(0),maxx-r);hi=min(a1,minx+r)
    if hi<=lo:continue
    ylo=max(F(0),(k+l2-1)*s-r);yhi=min(2*s,k*s+r)
    if yhi<=ylo:continue
    box=[lo,hi,ylo,yhi]
    for l in ls[2:]:
     if l%2:x,y=F(0),min(s/2,r-(l-1)*s/2)
     else:x,y=max(F(0),l*s/2-r),s/2
     if y<=x:break
     box.extend((x,y))
    else:raw.append(((j,k,ls),tuple(box)))
 unique={b:tag for tag,b in raw};kept=[]
 for b,tag in unique.items():
  budget()
  if not any(o!=b and all(o[2*q]<=b[2*q] and b[2*q+1]<=o[2*q+1] for q in range(5)) for o in unique):kept.append((tag,b))
 # Independent direct atom guard: enumerate the selected distinct point set
 # and verify both receiver endpoints against every point in every axis.
 atom_checks=0;max_abs_period_indices=[0]*5
 for (j,k,ls),b in kept:
  l1,l2=ls[:2];r=root(prod(ls),u,v)/2
  tails=[range(-(l-1)//2,(l-1)//2+1) if l%2 else range(1-l//2,1+l//2) for l in ls[2:]]
  atoms=[]
  for jj,kk,m3,m4,m5 in product(range(j,j+l1),range(k,k+l2),*tails):
   atoms.append((jj*a1+(h if kk%2 else F(0)),kk*s,m3*s,m4*s,m5*s))
   inds=(jj,kk//2,m3,m4,m5)
   max_abs_period_indices=[max(a,abs(z)) for a,z in zip(max_abs_period_indices,inds)]
  assert len(atoms)==prod(ls) and len(set(atoms))==prod(ls)
  for atom in atoms:
   for q,p in enumerate(atom):
    assert p-r<=b[2*q] and b[2*q+1]<=p+r, (j,k,ls,q)
    atom_checks+=1
  budget()
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
 for q in range(5):
  cuts=sorted({b[2*q+k] for b in retained for k in (0,1)});axis=[]
  for lo,hi in zip(cuts,cuts[1:]):
   mask=sum(1<<i for i,b in enumerate(retained) if b[2*q]<=lo and hi<=b[2*q+1])
   if mask:axis.append((hi-lo,mask))
  segments.append(axis)
 @lru_cache(None)
 def grid(q,mask):
  budget()
  if q==5:return F(1)
  out=F(0)
  for d,nxt in segments[q]:
   both=nxt if mask<0 else mask&nxt
   if both:out+=d*grid(q+1,both)
  return out
 assert volume==grid(0,-1);ratio=4*volume
 c=json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text());C=decode(F,c['coefficients']);diff=ratio-C
 lo,hi=bounds(ratio.a);dl,dh=bounds(diff.a)
 return {'status':'PASS_exact_improvement' if dl>0 else 'NO_IMPROVEMENT_exact','tau':str(tau),'row_shift':'tau*a1','input_blocks':1024,'start_ranges':{'j':[-6,6],'k':[-6,6]},
  'raw_boxes':len(raw),'retained_boxes':len(kept),'receiver_domain':'[0,a1] x [0,2s] x [0,s/2]^3','period_mass':2,'ratio_normalization':'4*absolute_union_volume',
  'witness_rule':'each retained box forces inclusion of N selected unit atoms at r=N^(1/5)/2',
  'direct_atom_endpoint_checks':atom_checks,'distinct_N_guard':True,'max_abs_selected_period_indices':max_abs_period_indices,
  'fix':'row interval [(k+l2-1)*s-r,k*s+r]','union_engines_agree':True,'coordinate_segment_counts':[len(x) for x in segments],
  'recursive_nodes':nodes,'wall_seconds':time.monotonic()-started,'ratio_coefficients':{f'{i//5},{i%5}':str(x) for i,x in enumerate(ratio.a) if x},
  'ratio_interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)],'difference_from_exact_C_interval':[str(dl),str(dh)],'strictly_above_exact_C':dl>0,
  'C_source_sha256':hashlib.sha256((ROOT/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest(),'limits':{'wall_seconds':LIMIT_SECONDS,'recursive_nodes':LIMIT_NODES},'field_source_sha256':hashlib.sha256((ROOT/'lower/certify_rectangular.py').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
def main():
 if len(sys.argv)!=2:raise SystemExit('usage: two_row_staggered.py tau')
 tau=Q(sys.argv[1]);assert tau in (Q(1,8),Q(1,4),Q(1,2))
 try:out=run(tau)
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','tau':str(tau),'reason':str(e)}
 (HERE/f'corrected_tau_{tau.numerator}_{tau.denominator}_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
