#!/usr/bin/env python3
"""Body-centered two-orbit finite witnesses, exact strict rational radii."""
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
 'anisotropic_9_8':[Q(7,4)*Q(9,8)]+[Q(17,16)*Q(9,8)]*4,
 'anisotropic_19_16':[Q(7,4)*Q(19,16)]+[Q(17,16)*Q(19,16)]*4,
 'anisotropic_5_4':[Q(7,4)*Q(5,4)]+[Q(17,16)*Q(5,4)]*4,
 'isotropic_11_8':[Q(11,8)]*5,
 'isotropic_23_16':[Q(23,16)]*5,
 'isotropic_3_2':[Q(3,2)]*5}
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

def indices(orbit,l):
 if orbit==0:return tuple(range(-(l-1)//2,(l-1)//2+1)) if l%2 else tuple(range(1-l//2,1+l//2))
 return tuple(range(-(l-1)//2,(l-1)//2+1)) if l%2 else tuple(range(-l//2,l//2))

def run(name):
 started=time.monotonic();nodes=0;configs=0
 def budget(node=False):
  nonlocal nodes
  if node:nodes+=1
  if nodes>NODE_CAP:raise RuntimeError('node cap')
  if time.monotonic()-started>SECONDS:raise TimeoutError('wall cap')
 a=CANDIDATES[name];domain=tuple(x/2 for x in a)
 # Coordinates are stored in units a_i/2. Orbit parity is all-even or all-odd.
 extrema={}
 for orbit in (0,1):
  for l in range(1,5):
   inds=indices(orbit,l);assert len(inds)==l
   extrema[orbit,l]=(2*min(inds)+orbit,2*max(inds)+orbit)
 raw={}
 def add(lengths0,lengths1):
  nonlocal configs
  configs+=1
  if configs%64==0:budget()
  n0=prod(lengths0) if lengths0 else 0;n1=prod(lengths1) if lengths1 else 0;n=n0+n1
  r=lower_side(n)/2;box=[]
  for q in range(5):
   lows=[];highs=[]
   for orbit,ls in ((0,lengths0),(1,lengths1)):
    if ls:
     low,high=extrema[orbit,ls[q]];lows.append(low);highs.append(high)
   lo=max(Q(0),max(highs)*a[q]/2-r);hi=min(domain[q],min(lows)*a[q]/2+r)
   if hi<=lo:return
   box.extend((lo,hi))
  raw[tuple(box)]=(lengths0,lengths1,n)
 for ls in product(range(1,5),repeat=5):add(ls,None);add(None,ls)
 pairs=tuple((i,j) for i in range(1,5) for j in range(1,5) if abs(i-j)<=1)
 for pairvec in product(pairs,repeat=5):add(tuple(p[0] for p in pairvec),tuple(p[1] for p in pairvec))
 budget();kept=[];all_boxes=list(raw)
 for b in all_boxes:
  budget()
  if not any(o!=b and all(o[2*q]<=b[2*q] and b[2*q+1]<=o[2*q+1] for q in range(5)) for o in all_boxes):kept.append(b)
 # Independent direct point guard for the retained witnesses.
 guards=0;max_abs_indices=[0]*5;max_selected=0
 for b in kept:
  l0,l1,n=raw[b];r=lower_side(n)/2;atoms=[]
  for orbit,ls in ((0,l0),(1,l1)):
   if ls:
    for inds in product(*(indices(orbit,l) for l in ls)):
     atoms.append(tuple(2*z+orbit for z in inds))
     max_abs_indices=[max(old,abs(z)) for old,z in zip(max_abs_indices,inds)]
  assert len(atoms)==n==len(set(atoms));max_selected=max(max_selected,n)
  assert (2*r)**5<n
  # Check extrema from the actual enumerated points, not the construction table.
  for q in range(5):
   lo=min(pt[q] for pt in atoms)*a[q]/2;hi=max(pt[q] for pt in atoms)*a[q]/2
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
 assert volume==grid(0,-1);ratio=16*volume
 cp=ROOT/'lower_round003/a2_equals_a3_results.json';C=json.loads(cp.read_text());cl,ch=map(Q,C['ratio_interval'])
 lower_improvement=ratio>ch;upper_failure=ratio<cl
 with localcontext() as ctx:
  ctx.prec=50;display=str(Decimal(ratio.numerator)/Decimal(ratio.denominator))
 return {'status':'PASS_exact_strict_candidate_above_C_pending_audit' if lower_improvement else ('NO_IMPROVEMENT_exact' if upper_failure else 'UNRESOLVED_C_interval'),
  'candidate':name,'dimension':5,'spacings':[str(x) for x in a],'period_mass':2,'orbits':'all coordinates even half-indices OR all coordinates odd half-indices','receiver_domain':'product [0,a_i/2]','ratio_normalization':'16*absolute_union_volume','finite_configurations_enumerated':configs,'distinct_raw_boxes':len(raw),'retained_boxes':len(kept),'union_engines_agree':True,'coordinate_segment_counts':[len(s) for s in segments],'distinct_selected_atom_guard':True,'direct_geometry_endpoint_inequalities':guards,'max_atoms_in_retained_witness':max_selected,'max_abs_selected_orbit_lattice_indices':max_abs_indices,'strict_radius_rule':'side=floor_lower_2^-80(N^(1/5)); side^5<N','ratio_exact_rational':str(ratio),'ratio_decimal_display_only':display,'strictly_above_C_upper':lower_improvement,'strictly_below_C_lower':upper_failure,'frozen_reference':'complete algebraic C_cert, rational isolating bounds used only in stated directions','C_source_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),'nodes':nodes,'wall_seconds':time.monotonic()-started,'limits':{'wall_seconds':SECONDS,'nodes':NODE_CAP},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
def main():
 if len(sys.argv)!=2 or sys.argv[1] not in CANDIDATES:raise SystemExit('choose '+','.join(CANDIDATES))
 name=sys.argv[1]
 try:out=run(name)
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','candidate':name,'reason':str(e),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/(name+'_results.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
