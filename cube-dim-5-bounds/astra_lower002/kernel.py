#!/usr/bin/env python3
"""Exact exposed-facet insertion response, with a registered one-axis cell."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import json,sys,time,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
F,u,v,bounds,_,_=field()
A=(1+u**3/2,)+((1+u)/2,)*4
ZERO=F(0);ONE=F(1)
def encode(x):return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.a) if c}
def decode(o):
 a=[Q(0)]*25
 for key,value in o.items():i,j=map(int,key.split(','));a[5*i+j]=Q(value)
 return F(a=a)
def disp(x):return [float(t) for t in bounds(x.a)]
@lru_cache(None)
def root(n):
 p=n;i=j=0
 while p%2==0:p//=2;i+=1
 while p%3==0:p//=3;j+=1
 assert p==1
 return u**i*v**j
def baseline():
 raw=[]
 for ls in product(range(1,5),repeat=5):
  n=prod(ls);r=root(n)/2;box=[]
  for aa,l in zip(A,ls):
   start=-(l-1)//2 if l%2 else 1-l//2
   lo=max(ZERO,(start+l-1)*aa-r);hi=min(aa/2,start*aa+r)
   if hi<lo:break
   box.extend((lo,hi))
  else:raw.append({'lengths':ls,'n':n,'r':r,'box':tuple(box)})
 positive=[row['box'] for row in raw if all(row['box'][2*j]<row['box'][2*j+1] for j in range(5))]
 keep=tuple(sorted(set(b for b in positive if not any(o!=b and all(o[2*j]<=b[2*j] and b[2*j+1]<=o[2*j+1] for j in range(5)) for o in positive))))
 return raw,keep
RAW,BASE=baseline()
def groups():
 out={}
 for row in RAW:
  b=row['box']
  for axis in range(5):
   if any(b[2*j]==b[2*j+1] for j in range(5) if j!=axis):continue
   for direction in (-1,1):
    t=b[2*axis+(direction==1)]
    if direction==-1 and t==ZERO:continue
    if direction==1 and t==A[axis]/2:continue
    key=(axis,t,direction)
    project=tuple(x for j in range(5) if j!=axis for x in b[2*j:2*j+2])
    out.setdefault(key,[]).append((project,row['n'],row['r']))
 result=[]
 for (axis,t,direction),active in out.items():
  block=[]
  for b in BASE:
   crossed=(b[2*axis]<=t and t<b[2*axis+1]) if direction==1 else (b[2*axis]<t and t<=b[2*axis+1])
   if crossed:block.append(tuple(x for j in range(5) if j!=axis for x in b[2*j:2*j+2]))
  # Remove entire candidate projections covered by a baseline blocker.
  active=[p for p in active if not any(all(o[2*j]<=p[0][2*j] and p[0][2*j+1]<=o[2*j+1] for j in range(4)) for o in block)]
  if active:result.append((axis,t,direction,active,tuple(block)))
 return result
GROUPS=groups()

def register_cell():
 target=A[1]/4;half=A[1]/2
 # Use all raw max/min alternatives, even from infeasible blocks. They safely
 # refine, rather than undercount, the displacement arrangement.
 radii=set();endpoints={ZERO,half}
 for ls in product(range(1,5),repeat=5):
  r=root(prod(ls))/2;radii.add(r);l=ls[1];start=-(l-1)//2 if l%2 else 1-l//2
  endpoints.update(((start+l-1)*A[1]-r,start*A[1]+r))
 events={ZERO,half}
 for e in endpoints:
  if ZERO<e and e<half:events.add(e)
  for r in radii:
   for k in range(-4,5):
    for value in (e+r-k*A[1],e-r-k*A[1]):
     for z in (value,-value):
      if ZERO<z and z<half:events.add(z)
 cuts=sorted(events)
 assert target not in events,'registered target lies exactly on event; stop for new registration'
 lower=max(e for e in cuts if e<target);upper=min(e for e in cuts if target<e)
 assert ZERO<lower and upper<F(Q(1,2))
 rows=[]
 for theta in (Q(1,4),Q(1,2),Q(3,4)):
  t=lower+(upper-lower)*F(theta)
  rows.append({'theta':str(theta),'position_coefficients':encode(t),'position_display_only':disp(t)})
 out={'status':'event_cell_registered_before_kernel_evaluation','coordinate':1,'lower_coefficients':encode(lower),'upper_coefficients':encode(upper),'lower_display_only':disp(lower),'upper_display_only':disp(upper),'target_display_only':disp(target),'event_count':len(events),'radii_count':len(radii),'raw_endpoint_alternatives':len(endpoints),'raw_feasible_or_degenerate_boxes':len(RAW),'positive_baseline_boxes':len(BASE),'exposed_face_groups':len(GROUPS),'isolated_satellite_halo_fraction':0,'isolated_halo_reason':'all selected t lie strictly in(0,1/2), remaining physical coordinates0, hence satellite lies in interior of the primary unit-atom receiver cube','nodes':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'event_cell.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

def kernel_half(z):
 begin=time.monotonic();nodes=0;patches=0;face_rows=[]
 def cap():
  if time.monotonic()-begin>120 or nodes>200000:raise TimeoutError('registered cap')
 @lru_cache(None)
 def count_at(axis,x,r):
  return sum(x-r<=k*A[axis]+z[axis] and k*A[axis]+z[axis]<=x+r for k in range(-4,5))
 @lru_cache(None)
 def split_counts(axis,lo,hi,r):
  cuts={lo,hi}
  for k in range(-4,5):
   atom=k*A[axis]+z[axis]
   for e in (atom-r,atom+r):
    if lo<e and e<hi:cuts.add(e)
  cuts=sorted(cuts);out=[]
  for left,right in zip(cuts,cuts[1:]):
   n=count_at(axis,(left+right)/2,r)
   if n:out.append((left,right,n))
  return tuple(out)
 @lru_cache(None)
 def integrate(positive,blockers):
  nonlocal nodes
  nodes+=1;cap()
  if not positive:return ZERO
  if blockers:
   positive=tuple(p for p in positive if not any(all(b[2*j]<=p[0][2*j] and p[0][2*j+1]<=b[2*j+1] for j in range(len(b)//2)) for b in blockers))
   if not positive:return ZERO
  cuts=sorted({x for b,w in positive for x in b[:2]}|{x for b in blockers for x in b[:2]})
  total=ZERO
  for lo,hi in zip(cuts,cuts[1:]):
   pp=tuple(sorted(set((b[2:],w) for b,w in positive if b[0]<=lo and hi<=b[1])))
   if not pp:continue
   bb=tuple(sorted(set(b[2:] for b in blockers if b[0]<=lo and hi<=b[1])))
   if not pp[0][0]:
    if not bb:total+=(hi-lo)*max(w for b,w in pp)
   else:total+=(hi-lo)*integrate(pp,bb)
  return total
 for axis,t,direction,active,blockers in GROUPS:
  cap();weighted=[];others=[j for j in range(5) if j!=axis]
  for project,n,r in active:
   face_count=count_at(axis,t,r)
   if not face_count:continue
   splits=[split_counts(j,project[2*d],project[2*d+1],r) for d,j in enumerate(others)]
   for item in product(*splits):
    count=face_count*prod(x[2] for x in item)
    weight=r*count/(5*n)
    box=tuple(y for x in item for y in x[:2]);weighted.append((box,weight))
  weighted=tuple(sorted(set(weighted)));patches+=len(weighted)
  value=integrate(weighted,blockers)
  if value!=ZERO:face_rows.append({'axis':axis,'direction':direction,'plane':encode(t),'contribution':encode(value),'display_only':disp(value)})
 total=sum((decode(row['contribution']) for row in face_rows),ZERO)
 return {'status':'exact_halfcell_exposed_kernel','shift_coefficients':[encode(x) for x in z],'kernel_coefficients':encode(total),'kernel_display_only':disp(total),'axis_contributions':[encode(sum((decode(row['contribution']) for row in face_rows if row['axis']==j),ZERO)) for j in range(5)],'nonzero_faces':face_rows,'weighted_patches':patches,'recursive_nodes':nodes,'seconds':time.monotonic()-begin,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
def evaluate(index,sign):
 reg=json.loads((HERE/'event_cell.json').read_text());t=decode(reg['nodes'][index]['position_coefficients']);z=[ZERO]*5;z[1]=t*sign
 try:out=kernel_half(tuple(z))
 except TimeoutError as e:out={'status':'COMPUTE_LIMIT','reason':str(e)}
 out.update(node=index,reflection_sign=sign)
 (HERE/f'kernel_{index}_{"plus" if sign>0 else "minus"}.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ('nonzero_faces','shift_coefficients','kernel_coefficients','axis_contributions')},indent=2))
def finish():
 reg=json.loads((HERE/'event_cell.json').read_text());C=decode(json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients']);values=[]
 for i in range(3):
  plus=json.loads((HERE/f'kernel_{i}_plus.json').read_text());minus=json.loads((HERE/f'kernel_{i}_minus.json').read_text());assert plus['status']==minus['status']=='exact_halfcell_exposed_kernel'
  B=16*(decode(plus['kernel_coefficients'])+decode(minus['kernel_coefficients']));D=B-C
  values.append(D)
 assert 2*values[1]==values[0]+values[2]
 left=2*values[0]-values[1];right=2*values[2]-values[1]
 def entry(x):return {'coefficients':encode(x),'interval':[str(q) for q in bounds(x.a)],'display_only':disp(x)}
 out={'status':'STRICT_NEGATIVE_ON_REGISTERED_OPEN_CELL' if left<ZERO and right<ZERO else 'OTHER_SIGN_REQUIRES_REVIEW','normalization':'D(z)=full numerator boundary derivative-C; isolated halo contribution0 on this cell','node_derivatives':[entry(x) for x in values],'left_limit':entry(left),'right_limit':entry(right),'affine_midpoint_identity_exact':True,'cell_scope':'open physical displacement interval in event_cell.json, z=(0,t,0,0,0); not full5D displacement exclusion','C_coefficients_sha256':hashlib.sha256((ROOT/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest()}
 (HERE/'cell_derivative_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('node_derivatives','left_limit','right_limit')},indent=2));print('endpoint limits:',disp(left),disp(right))
if __name__=='__main__':
 if sys.argv[1]=='register':register_cell()
 elif sys.argv[1]=='evaluate':evaluate(int(sys.argv[2]),int(sys.argv[3]))
 elif sys.argv[1]=='finish':finish()
