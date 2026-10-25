"""Independent inward-rational witness boxes; integer union, no number-field code."""
from fractions import Fraction as Q
from itertools import product
from functools import lru_cache
from pathlib import Path
import json,math,time,hashlib
ROOT=Path("/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds")
D=10**50; GRID=10**35

def rootfloor(m):
 lo=0;hi=10*D;target=m*D**6
 while hi-lo>1:
  mid=(lo+hi)//2
  if mid**6<=target:lo=mid
  else:hi=mid
 assert lo**6<=target<(lo+1)**6
 return Q(lo,D),Q(lo+1,D)

def ceilq(x):return -((-x.numerator)//x.denominator)

@lru_cache(None)
def union(B):
 if not B:return 0
 if len(B[0])==2:
  pairs=sorted(B);start,end=pairs[0];vol=0
  for left,right in pairs[1:]:
   if left>end:vol+=end-start;start,end=left,right
   else:end=max(end,right)
  return vol+end-start
 cuts=sorted({x for b in B for x in b[:2]});vol=0
 for left,right in zip(cuts,cuts[1:]):
  child=tuple(sorted({b[2:] for b in B if b[0]<=left and right<=b[1]}))
  vol+=(right-left)*union(child)
 return vol

start=time.time();u0,u1=rootfloor(2);v0,v1=rootfloor(3)
spacing=[(1+u0**4/2,1+u1**4/2),((2+u0**5*v0)/4,(2+u1**5*v1)/4)]+[((1+u0)/2,(1+u1)/2)]*4
spacing[1]=(u0**2*v0**2/3+u0**2*v0**4/6,u1**2*v1**2/3+u1**2*v1**4/6)
boxes=[];witness=[];margin=Q(1)
for ell in product(range(1,5),repeat=6):
 mass=math.prod(ell);s0,s1=rootfloor(mass);rad=s0/2;box=[]
 for l,(al,ah) in zip(ell,spacing):
  if l%2:low=Q(0);high=min(al/2,rad-(l-1)*ah/2)
  else:low=max(Q(0),l*ah/2-rad);high=al/2
  # Strict inward rounding: positive margin even for integral grid endpoints.
  lowint=ceilq(low*GRID)+1;highint=(high*GRID).numerator//(high*GRID).denominator-1
  if highint<=lowint:break
  box.extend((lowint,highint))
 else:
  for j,(l,(al,ah)) in enumerate(zip(ell,spacing)):
   first=-(l//2) if l%2 else 1-l//2
   for index in (first,first+l-1):
    atomlo,atomhi=sorted((index*al,index*ah))
    assert Q(box[2*j],GRID)-atomhi>-rad
    assert Q(box[2*j+1],GRID)-atomlo<rad
    margin=min(margin,rad+Q(box[2*j],GRID)-atomhi,rad-Q(box[2*j+1],GRID)+atomlo)
   assert 0<Q(box[2*j],GRID)<Q(box[2*j+1],GRID)<al/2
  assert (2*rad)**6<=mass
  boxes.append(tuple(box));witness.append(ell)
unique=set(boxes)
kept=tuple(sorted(b for b in unique if not any(c!=b and all(c[2*j]<=b[2*j] and b[2*j+1]<=c[2*j+1] for j in range(6)) for c in unique)))
R=Q(64*union(kept),GRID**6)
target=Q(2067397333267,10**12)
assert R>target
threshold=(R+target)/(2*R);K=1
while threshold*R*Q(2*K+1,2*K+5)**6<=target:K*=2
finite_ratio=threshold*R*Q(2*K+1,2*K+5)**6
assert margin>0 and finite_ratio>target
opening=json.loads((ROOT/'rounds/001/n5_readonly_snapshot.json').read_text())
assert R>Q(opening['interval'][1])
old=json.loads((ROOT/'results/bounds_ledger.json').read_text())
frozen=json.loads((ROOT/'rounds/003/frozen_lower.json').read_text())
old_upper=Q(frozen['lower']['exact_interval'][1])
assert R>old_upper and finite_ratio>old_upper
current=json.loads((ROOT/'rounds/008/frozen_lower_before_candidate.json').read_text())
assert R>Q(current['lower']['exact_interval'][1])
assert finite_ratio>Q(current['lower']['exact_interval'][1])
r4=old['inherited_lower']['interval'][1]
assert R>Q(r4)
result={'status':'pass','method':'strict inward rational boxes and integer union; no author imports or algebraic field arithmetic','grid_denominator':GRID,'root_denominator':D,'candidate_boxes':len(boxes),'retained_boxes':len(kept),'exact_rational_lower':str(R),'simple_strict_lower':str(target),'spacing_change':'a2=u^2*v^2/3+u^2*v^4/6','strictly_above_full_previous_R6_upper':True,'old_R6_upper':str(old_upper),'strictly_above_current_full_upper':True,'current_full_upper':current['lower']['exact_interval'][1],'finite_L1_witness':{'source_index_radius':K+2,'receiver_index_radius':K,'threshold':str(threshold),'bump_halfwidth':str(margin/2),'bump':'uniform density (2 eta)^-6 on [-eta,eta]^6','finite_ratio_strictly_above_target':bool(finite_ratio>target)},'strictly_above_full_R5_upper':True,'R5_upper':opening['interval'][1],'strictly_above_full_R4_upper':True,'R4_upper':r4,'elapsed_seconds':time.time()-start,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('inward_smallstep_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('R5_upper','R4_upper')},indent=2))
