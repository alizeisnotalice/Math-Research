#!/usr/bin/env python3
from pathlib import Path
from itertools import product,combinations_with_replacement
from fractions import Fraction as Q
import json,sys,hashlib,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'astra_lower002'))
import kernel as K
F=K.F;A=K.A;ZERO=K.ZERO
begin=time.monotonic()
def arrangement(axis):
 cuts=sorted({x for b in K.BASE for x in b[2*axis:2*axis+2]}|{ZERO,A[axis]/2})
 return [(lo,hi,sum(1<<j for j,b in enumerate(K.BASE) if b[2*axis]<=lo and hi<=b[2*axis+1])) for lo,hi in zip(cuts,cuts[1:])]
first=arrangement(0);tail=arrangement(1);tail_axes=[arrangement(j) for j in range(1,5)]
assert all([(lo,hi) for lo,hi,m in rows]==[(lo,hi) for lo,hi,m in tail] for rows in tail_axes)
multis=list(combinations_with_replacement(range(len(tail)),4));multis.sort(key=lambda t:(len(set(t)),t))
selected=None;checked=0
for tails in multis:
 for fi in range(len(first)):
  checked+=1;mask=first[fi][2]
  for axis,j in enumerate(tails):mask&=tail_axes[axis][j][2]
  if not mask:
   selected=(fi,tails);break
 if selected:break
assert selected
coarse=[first[selected[0]]]+[tail_axes[axis][j] for axis,j in enumerate(selected[1])]
def event_cell(axis,lo,hi):
 target=(lo+hi)/2;half=A[axis]/2;radii={row['r'] for row in K.RAW};endpoints={ZERO,half}
 for lengths in product(range(1,5),repeat=5):
  r=K.root(__import__('math').prod(lengths))/2;l=lengths[axis];start=-(l-1)//2 if l%2 else 1-l//2
  endpoints.update(((start+l-1)*A[axis]-r,start*A[axis]+r));radii.add(r)
 cuts={lo,hi}
 for e in endpoints:
  if lo<e and e<hi:cuts.add(e)
  for r in radii:
   for integer in range(-4,5):
    for t in (e+r-integer*A[axis],e-r-integer*A[axis]):
     for value in (t,-t):
      if lo<value and value<hi:cuts.add(value)
 cuts=sorted(cuts)
 if target in cuts:
  i=cuts.index(target)
  target=(cuts[i]+cuts[i+1])/2 if i+1<len(cuts) else (cuts[i-1]+cuts[i])/2
 left=max(x for x in cuts if x<target);right=min(x for x in cuts if target<x)
 return left,right,(left+right)/2,len(cuts)
cache={};fine=[]
for axis,(lo,hi,mask) in enumerate(coarse):
 key=(0 if axis==0 else 1,lo,hi)
 if key not in cache:cache[key]=event_cell(axis,lo,hi)
 fine.append(cache[key])
z=tuple(row[2] for row in fine)
assert not any(all(b[2*j]<=z[j] and z[j]<=b[2*j+1] for j in range(5)) for b in K.BASE)
assert all(ZERO<row[0] and row[1]<A[j]/2 for j,row in enumerate(fine))
zero=[row['lengths'] for row in K.RAW if any(row['box'][2*j]==row['box'][2*j+1] for j in range(5))]
facets={}
for row in K.RAW:
 for axis in range(5):
  for direction in (-1,1):
   t=row['box'][2*axis+(direction==1)]
   if t==ZERO or t==A[axis]/2:continue
   facets.setdefault((axis,t,direction),set()).add(row['r'])
unequal=[key for key,radii in facets.items() if len(radii)>1]
out={'status':'SELECTED_UNCOVERED_OPEN_CELL_AND_PASSED_CONCAVITY_GEOMETRY' if not zero and not unequal else 'STRUCTURAL_AUDIT_FAILED_STOP','selection':{'coarse_first_index':selected[0],'coarse_tail_indices':selected[1],'tested_mask_patterns':checked,'tail_distinct_intervals':len(set(selected[1]))},'coarse_box':[[K.encode(lo),K.encode(hi)] for lo,hi,mask in coarse],'fine_cell':[[K.encode(lo),K.encode(hi)] for lo,hi,t,n in fine],'shift_coefficients':[K.encode(t) for t in z],'shift_display_only':[K.disp(t) for t in z],'event_counts_within_coarse_axes':[row[3] for row in fine],'kappa':1,'cell_exactly_disjoint_from_baseline_union':True,'zero_width_primary_boxes':zero,'unequal_radius_same_oriented_plane_groups':len(unequal),'same_oriented_plane_groups':len(facets),'raw_boxes':len(K.RAW),'seconds':time.monotonic()-begin,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'kernel_scalar_geometry_source_sha256':hashlib.sha256((HERE.parent/'astra_lower002/kernel.py').read_bytes()).hexdigest()}
(HERE/'selected_cell.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('coarse_box','fine_cell','shift_coefficients')},indent=2))
