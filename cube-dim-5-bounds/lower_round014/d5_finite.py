#!/usr/bin/env python3
"""D5 full-five-axis parity lattice witnesses, exact strict rational radii."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
DEN=2**80; SECONDS=120; NODE_CAP=200000
CANDIDATES={
 'anisotropic_27_32':[Q(7,4)*Q(27,32)]+[Q(17,16)*Q(27,32)]*4,
 'anisotropic_7_8':[Q(7,4)*Q(7,8)]+[Q(17,16)*Q(7,8)]*4,
 'anisotropic_29_32':[Q(7,4)*Q(29,32)]+[Q(17,16)*Q(29,32)]*4,
 'isotropic_31_32':[Q(31,32)]*5,
 'isotropic_1':[Q(1)]*5,
 'isotropic_33_32':[Q(33,32)]*5}
@lru_cache(None)
def lower_side(n):
 lo,hi=0,8*DEN
 while lo+1<hi:
  mid=(lo+hi)//2
  if mid**5<=n*DEN**5:lo=mid
  else:hi=mid
 if lo**5==n*DEN**5:lo-=1
 out=Q(lo,DEN);assert out**5<n and (out+Q(1,DEN))**5>=n
 return out

def run(name):
 started=time.monotonic();nodes=0;configs=0
 def budget(node=False):
  nonlocal nodes
  if node:nodes+=1
  if nodes>NODE_CAP:raise RuntimeError('node cap')
  if time.monotonic()-started>SECONDS:raise TimeoutError('wall cap')
 a=CANDIDATES[name];domain=tuple(a)
 # Six nearest-block choices per axis on [0,a_i]. This is a registered
 # finite family rather than an assertion that every relevant block is used.
 choices=((1,0),(1,1),(2,0),(3,-1),(3,0),(4,-1))
 raw={}
 for choices5 in product(choices,repeat=5):
  configs+=1
  if configs%32==0:budget()
  ls=tuple(x[0] for x in choices5);ks=tuple(x[1] for x in choices5)
  allodd=int(all(l%2 for l in ls));sign=1 if sum(ks)%2==0 else -1
  n=(prod(ls)+sign*allodd)//2
  if not n:continue
  r=lower_side(n)/2;box=[]
  for q in range(5):
   if all(ls[t]==1 for t in range(5) if t!=q):
    low=ks[q]+(sum(ks)%2)
    high=ks[q]+ls[q]-1-((sum(ks)+ls[q]-1)%2)
   else:low,high=ks[q],ks[q]+ls[q]-1
   assert low<=high
   lo=max(Q(0),high*a[q]-r);hi=min(domain[q],low*a[q]+r)
   if hi<=lo:break
   box.extend((lo,hi))
  else:raw[tuple(box)]=(ls,ks,n)
 budget();kept=[];all_boxes=list(raw)
 for b in all_boxes:
  budget()
  if not any(o!=b and all(o[2*q]<=b[2*q] and b[2*q+1]<=o[2*q+1] for q in range(5)) for o in all_boxes):kept.append(b)
 guards=0;max_abs_indices=[0]*5;max_selected=0
 for b in kept:
  ls,ks,n=raw[b];r=lower_side(n)/2
  atoms=[tuple(pt) for pt in product(*(range(k,k+l) for k,l in zip(ks,ls))) if sum(pt)%2==0]
  assert len(atoms)==n==len(set(atoms));max_selected=max(max_selected,n)
  assert (2*r)**5<n
  for atom in atoms:max_abs_indices=[max(old,abs(z)) for old,z in zip(max_abs_indices,atom)]
  for q in range(5):
   lo=min(pt[q] for pt in atoms)*a[q];hi=max(pt[q] for pt in atoms)*a[q]
   assert hi-r<=b[2*q] and b[2*q+1]<=lo+r
   guards+=2
  budget()
 retained=tuple(sorted(kept))
 @lru_cache(None)
 def sweep(boxes):
  budget(True)
  if not boxes:return Q(0)
  if len(boxes[0])==2:
   out=Q(0);lo=hi=None
   for x,y in sorted(boxes):
    if lo is None:lo,hi=x,y
    elif hi<x:out+=hi-lo;lo,hi=x,y
    else:hi=max(hi,y)
   return out+hi-lo
  cuts=sorted({x for b in boxes for x in b[:2]});out=Q(0)
  for lo,hi in zip(cuts,cuts[1:]):
   out+=(hi-lo)*sweep(tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]})))
  return out
 volume=sweep(retained);segments=[]
 for q in range(5):
  budget();cuts=sorted({b[2*q+k] for b in retained for k in (0,1)});axis=[]
  for lo,hi in zip(cuts,cuts[1:]):
   mask=sum(1<<i for i,b in enumerate(retained) if b[2*q]<=lo and hi<=b[2*q+1])
   if mask:axis.append((hi-lo,mask))
  segments.append(axis)
 @lru_cache(None)
 def grid(q,mask):
  budget(True)
  if q==5:return Q(1)
  out=Q(0)
  for d,other in segments[q]:
   nxt=other if mask<0 else mask&other
   if nxt:out+=d*grid(q+1,nxt)
  return out
 assert volume==grid(0,-1);ratio=2*volume
 cp=ROOT/'lower_round003/a2_equals_a3_results.json';C=json.loads(cp.read_text());cl,ch=map(Q,C['ratio_interval'])
 lower_improvement=ratio>ch;upper_failure=ratio<cl
 with localcontext() as ctx:
  ctx.prec=50;display=str(Decimal(ratio.numerator)/Decimal(ratio.denominator))
 return {'status':'PASS_exact_strict_candidate_above_C_pending_audit' if lower_improvement else ('NO_IMPROVEMENT_exact' if upper_failure else 'UNRESOLVED_C_interval'),
  'candidate':name,'dimension':5,'spacings':[str(x) for x in a],'period_mass':16,'orbits':'integer lattice coordinates j with sum(j) even; 16 residue classes mod2','receiver_domain':'product [0,a_i]','ratio_normalization':'2*absolute_union_volume','finite_configurations_enumerated':configs,'distinct_raw_boxes':len(raw),'retained_boxes':len(kept),'union_engines_agree':True,'coordinate_segment_counts':[len(s) for s in segments],'distinct_selected_atom_guard':True,'direct_geometry_endpoint_inequalities':guards,'max_atoms_in_retained_witness':max_selected,'max_abs_selected_integer_lattice_indices':max_abs_indices,'strict_radius_rule':'side=floor_lower_2^-80(N^(1/5)); side^5<N','ratio_exact_rational':str(ratio),'ratio_decimal_display_only':display,'strictly_above_C_upper':lower_improvement,'strictly_below_C_lower':upper_failure,'frozen_reference':'complete algebraic C_cert, rational isolating bounds used only in stated directions','C_source_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),'nodes':nodes,'wall_seconds':time.monotonic()-started,'limits':{'wall_seconds':SECONDS,'nodes':NODE_CAP},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
def main():
 if len(sys.argv)!=2 or sys.argv[1] not in CANDIDATES:raise SystemExit('choose '+','.join(CANDIDATES))
 name=sys.argv[1]
 try:out=run(name)
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','candidate':name,'reason':str(e),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/(name+'_results.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
