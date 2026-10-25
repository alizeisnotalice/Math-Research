#!/usr/bin/env python3
"""Finite nonperiodic two-atom cluster screen; exact rational enclosures."""
from fractions import Fraction as Q
from itertools import combinations
from functools import lru_cache
from pathlib import Path
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent
D=5; DEN=2**80; SECONDS=120; NODES=200000
CANDIDATES={
 'equal_axis':([Q(1,2),Q(1,2)],[Q(15,16),Q(0),Q(0),Q(0),Q(0)]),
 'equal_diagonal':([Q(1,2),Q(1,2)],[Q(1,2)]*5),
 'unequal_axis':([Q(31,32),Q(1,32)],[Q(15,16),Q(0),Q(0),Q(0),Q(0)])}
def root_bounds(m):
 lo,hi=0,2*DEN
 while lo+1<hi:
  mid=(lo+hi)//2
  if mid**D*m.denominator<=m.numerator*DEN**D:lo=mid
  else:hi=mid
 if lo**D*m.denominator==m.numerator*DEN**D:lo-=1;hi-=1
 a,b=Q(lo,DEN),Q(hi,DEN)
 assert a**D<m<=b**D
 return a,b

def run(name):
 started=time.monotonic();nodes=0
 def budget():
  nonlocal nodes
  nodes+=1
  if nodes>NODES:raise RuntimeError('node cap')
  if time.monotonic()-started>SECONDS:raise TimeoutError('wall cap')
 weights,d=CANDIDATES[name];points=[(Q(0),)*D,tuple(d)];assert len(set(points))==2
 assert all(w>0 for w in weights) and sum(weights)==1
 sides={};lower=[];upper=[];guards=0
 for size in (1,2):
  for subset in combinations(range(2),size):
   budget();mass=sum(weights[i] for i in subset);a,b=root_bounds(mass)
   sides[','.join(map(str,subset))]={'mass':str(mass),'strict_lower_side':str(a),'outer_side':str(b)}
   for side,boxes,strict in ((a,lower,True),(b,upper,False)):
    r=side/2;box=[]
    for q in range(D):
     vals=[points[i][q] for i in subset]
     box.extend((max(vals)-r,min(vals)+r))
    if all(box[2*q]<box[2*q+1] for q in range(D)):
     for i in subset:
      for q in range(D):
       assert points[i][q]-r<=box[2*q] and box[2*q+1]<=points[i][q]+r
       guards+=1
     if strict:assert (2*r)**D<mass
     boxes.append(tuple(box))
 def union(boxes):
  @lru_cache(None)
  def sweep(bs):
   budget()
   if not bs:return Q(0)
   if len(bs[0])==2:
    out=Q(0);lo=hi=None
    for x,y in sorted(bs):
     if lo is None:lo,hi=x,y
     elif hi<x:out+=hi-lo;lo,hi=x,y
     else:hi=max(hi,y)
    return out+hi-lo
   cuts=sorted({x for b in bs for x in b[:2]});out=Q(0)
   for lo,hi in zip(cuts,cuts[1:]):out+=(hi-lo)*sweep(tuple(sorted({b[2:] for b in bs if b[0]<=lo and hi<=b[1]})))
   return out
  answer=sweep(tuple(sorted(boxes)))
  # Algebraically distinct exact inclusion-exclusion for at most three boxes.
  second=Q(0)
  for k in range(1,len(boxes)+1):
   for inds in combinations(range(len(boxes)),k):
    budget();vol=Q(1)
    for q in range(D):vol*=max(Q(0),min(boxes[i][2*q+1] for i in inds)-max(boxes[i][2*q] for i in inds))
    second+=(-1)**(k+1)*vol
  assert answer==second
  return answer
 vlo,vhi=union(lower),union(upper)
 assert vlo<=vhi<Q(3,2)
 if name=='equal_axis':assert vlo<=Q(17,16)<=vhi
 if name=='equal_diagonal':assert vlo<=1<=vhi
 return {'status':'PASS_exact_full_level_enclosure_below_pair_obstruction','candidate':name,'dimension':D,'points':[[str(x) for x in p] for p in points],'weights':[str(w) for w in weights],'total_mass':'1','level':'1','subset_witnesses':sides,'union_boxes_lower_upper':[len(lower),len(upper)],'common_radius_rule':'one r=(subset_mass)^(1/5)/2 for all coordinates; rational inner/outer brackets','full_atomic_superlevel_measure_interval':[str(vlo),str(vhi)],'decimal_display_only':[float(vlo),float(vhi)],'strictly_below_3_over_2':True,'exact_pair_obstruction':'3/2','union_engines_agree':True,'direct_selected_atom_coordinate_checks':guards,'nodes':nodes,'wall_seconds':time.monotonic()-started,'limits':{'wall_seconds':SECONDS,'nodes':NODES},'root_denominator':str(DEN),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
def main():
 if len(sys.argv)!=2 or sys.argv[1] not in CANDIDATES:raise SystemExit('choose '+','.join(CANDIDATES))
 name=sys.argv[1]
 try:out=run(name)
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','candidate':name,'reason':str(e)}
 (HERE/(name+'_results.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
