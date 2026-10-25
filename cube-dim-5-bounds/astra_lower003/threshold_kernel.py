#!/usr/bin/env python3
"""Exact free-threshold hinge representation at one genuine halo-hole cell."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache,cmp_to_key
from itertools import product
from math import prod,comb
import json,hashlib,time,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'astra_lower002'))
import kernel as K
F=K.F;ZERO=K.ZERO;ONE=K.ONE;A=K.A
def cell():
 data=json.loads((HERE/'selected_cell.json').read_text());assert data['status']=='SELECTED_UNCOVERED_OPEN_CELL_AND_PASSED_CONCAVITY_GEOMETRY'
 return data,tuple(K.decode(o) for o in data['shift_coefficients'])
def projections(axis,b):return tuple(x for j in range(5) if j!=axis for x in b[2*j:2*j+2])
def prune(boxes):
 boxes=set(boxes)
 return tuple(sorted(b for b in boxes if not any(o!=b and all(o[2*j]<=b[2*j] and b[2*j+1]<=o[2*j+1] for j in range(len(b)//2)) for o in boxes)))
def contact_geometry():
 facets={}
 for row in K.RAW:
  for axis in range(5):
   for direction in (-1,1):
    t=row['box'][2*axis+(direction==1)]
    if t==ZERO or t==A[axis]/2:continue
    facets.setdefault((axis,t),{}).setdefault(direction,[]).append((projections(axis,row['box']),row['r'],row['n']))
 result=[]
 for (axis,t),dirs in facets.items():
  if -1 not in dirs or 1 not in dirs:continue
  left=dirs[1];right=dirs[-1]
  assert len({r for b,r,n in left})==len({r for b,r,n in right})==1
  rL,nL=left[0][1:];rR,nR=right[0][1:];intersections=[]
  for b,_,_ in left:
   for c,_,_ in right:
    q=tuple(x for j in range(4) for x in (max(b[2*j],c[2*j]),min(b[2*j+1],c[2*j+1])))
    if all(q[2*j]<q[2*j+1] for j in range(4)):intersections.append(q)
  blockers=prune(projections(axis,b) for b in K.BASE if b[2*axis]<t and t<b[2*axis+1])
  intersections=prune(b for b in intersections if not any(all(c[2*j]<=b[2*j] and b[2*j+1]<=c[2*j+1] for j in range(4)) for c in blockers))
  if intersections:result.append((axis,t,rL,nL,rR,nR,intersections,blockers))
 return result
CONTACTS=contact_geometry()
def volume_engine():
 begin=time.monotonic();counter=[0]
 @lru_cache(None)
 def volume(pos,neg=()):
  counter[0]+=1
  if counter[0]>200000 or time.monotonic()-begin>120:raise TimeoutError('registered volume cap')
  if not pos:return ZERO
  if neg:
   pos=tuple(p for p in pos if not any(all(n[2*j]<=p[2*j] and p[2*j+1]<=n[2*j+1] for j in range(len(n)//2)) for n in neg))
   if not pos:return ZERO
  cuts=sorted({x for b in pos for x in b[:2]}|{x for b in neg for x in b[:2]});total=ZERO
  for lo,hi in zip(cuts,cuts[1:]):
   pp=tuple(sorted({b[2:] for b in pos if b[0]<=lo and hi<=b[1]}))
   if not pp:continue
   nn=tuple(sorted({b[2:] for b in neg if b[0]<=lo and hi<=b[1]}))
   if not pp[0]:
    if not nn:total+=hi-lo
   else:total+=(hi-lo)*volume(pp,nn)
  return total
 return volume,counter,begin
def geometry():
 vol,count,begin=volume_engine();outward=ZERO;contact=ZERO;planes=[]
 for axis,t,direction,active,blockers in K.GROUPS:
  radii={r for b,n,r in active};assert len(radii)==1
  r=next(iter(radii));area=vol(prune(b for b,n,r in active),prune(blockers));outward+=area*r/5
 for axis,t,rL,nL,rR,nR,inter,block in CONTACTS:
  area=vol(inter,block)
  if area!=ZERO:
   contact+=area*(rL+rR)/5
   planes.append({'axis':axis,'plane':K.encode(t),'left_radius':K.encode(rL),'right_radius':K.encode(rR),'left_N':nL,'right_N':nR,'area':K.encode(area)})
 C=K.decode(json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients'])
 S=C-32*outward;total=32*contact
 assert ZERO<S and S<total
 out={'status':'FREE_THRESHOLD_CONCAVITY_GEOMETRY_CERTIFIED','outer_K_full':K.encode(32*outward),'initial_h_slope':K.encode(S),'total_contact_slope_drop':K.encode(total),'display_only':{'K':K.disp(32*outward),'S':K.disp(S),'contact_drop':K.disp(total)},'contact_planes':planes,'contact_plane_count':len(planes),'nodes':count[0],'seconds':time.monotonic()-begin,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'threshold_geometry.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('outer_K_full','initial_h_slope','total_contact_slope_drop','contact_planes')},indent=2))
def reflection_spec(index):
 specs=list(product(range(2),range(3),range(3)));f,p,q=specs[index]
 signs=(-1 if f else 1,)+(-1,)*p+(1,)*(2-p)+(-1,)*q+(1,)*(2-q)
 return signs,comb(2,p)*comb(2,q)
def run_reflection(index):
 original,z=cell();signs,multiplicity=reflection_spec(index);z=tuple(t*s for t,s in zip(z,signs))
 # Positive h=0 boundary intercept, reusing the already verified exposed
 # kernel. Threshold contact losses are separately reconstructed below.
 outer=K.kernel_half(z)
 vol,count,begin=volume_engine();hist=[];partition_checks=0
 @lru_cache(None)
 def at(axis,x,r):return sum(x-r<=integer*A[axis]+z[axis] and integer*A[axis]+z[axis]<=x+r for integer in range(-4,5))
 @lru_cache(None)
 def split(axis,lo,hi,rL,rR):
  cuts={lo,hi}
  for integer in range(-4,5):
   point=integer*A[axis]+z[axis]
   for r in (rL,rR):
    for e in (point-r,point+r):
     if lo<e and e<hi:cuts.add(e)
  cuts=sorted(cuts)
  return tuple((lo,hi,at(axis,(lo+hi)/2,rL),at(axis,(lo+hi)/2,rR)) for lo,hi in zip(cuts,cuts[1:]))
 for axis,t,rL,nL,rR,nR,intersections,blockers in CONTACTS:
  grouped={};other=[j for j in range(5) if j!=axis]
  baseL=at(axis,t,rL);baseR=at(axis,t,rR)
  for b in intersections:
   bands=[split(j,b[2*d],b[2*d+1],rL,rR) for d,j in enumerate(other)]
   for pieces in product(*bands):
    satelliteL=baseL*prod(p[2] for p in pieces);satelliteR=baseR*prod(p[3] for p in pieces)
    box=tuple(q for p in pieces for q in p[:2]);grouped.setdefault((satelliteL,satelliteR),[]).append(box)
  area_sum=ZERO
  for (sL,sR),boxes in grouped.items():
   area=vol(tuple(sorted(set(boxes))),blockers)
   if area==ZERO:continue
   alpha=(rL+rR)/5;beta=rL*sL/(5*nL)+rR*sR/(5*nR)
   hist.append({'axis':axis,'satellite_counts':[sL,sR],'primary_counts':[nL,nR],'area':K.encode(area),'alpha':K.encode(alpha),'beta':K.encode(beta)})
   area_sum+=area
  assert area_sum==vol(intersections,blockers);partition_checks+=1
 out={'status':'exact_reflected_hinge_histogram','reflection_index':index,'signs':signs,'multiplicity':multiplicity,'outer_intercept':outer['kernel_coefficients'],'outer_kernel_nodes':outer['recursive_nodes'],'outer_kernel_seconds':outer['seconds'],'hinges':hist,'contact_partition_checks':partition_checks,'contact_nodes':count[0],'contact_seconds':time.monotonic()-begin,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'outer_kernel_source_sha256':outer['source_sha256']}
 (HERE/f'reflection_{index:02}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('hinges','outer_intercept')},indent=2))
def finalize():
 info,z=cell();geom=json.loads((HERE/'threshold_geometry.json').read_text());S=K.decode(geom['initial_h_slope']);C=K.decode(json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients']);rows=[];B=ZERO;mult=0
 for index in range(18):
  data=json.loads((HERE/f'reflection_{index:02}.json').read_text());m=data['multiplicity'];mult+=m;B+=m*K.decode(data['outer_intercept'])
  for h in data['hinges']:
   area=m*K.decode(h['area']);alpha=K.decode(h['alpha']);beta=K.decode(h['beta']);rows.append((beta,alpha,area*alpha,area*beta))
 assert mult==32
 def compare(a,b):
  d=a[0]*b[1]-b[0]*a[1]
  return 0 if d==ZERO else -1 if d<ZERO else 1
 rows.sort(key=cmp_to_key(compare));events=[]
 for b,a,w,t in rows:
  if events and compare((b,a),(events[-1][0],events[-1][1]))==0:
   old=events[-1];events[-1]=(old[0],old[1],old[2]+w,old[3]+t)
  else:events.append((b,a,w,t))
 assert sum((e[2] for e in events),ZERO)==K.decode(geom['total_contact_slope_drop'])
 cumulative=ZERO;optimum=None
 for i,(b,a,w,t) in enumerate(events):
  if S<=cumulative+w:optimum=i;break
  cumulative+=w
 assert optimum is not None
 indices=sorted({max(0,optimum-1),optimum,min(len(events)-1,optimum+1)})
 D0=B+ONE-C
 def quotient_interval(num,den):
  nl,nh=K.bounds(num.a);dl,dh=K.bounds(den.a);assert dl>0
  choices=[nl/dl,nl/dh,nh/dl,nh/dh]
  return min(choices),max(choices)
 candidates=[]
 for i in indices:
  hnum,hden,_,_=events[i]
  numerator=D0*hden+S*hnum
  for b,a,w,t in events:
   hinge=w*hnum-t*hden
   if ZERO<hinge:numerator-=hinge
  hbounds=quotient_interval(hnum,hden);dbounds=quotient_interval(numerator,hden)
  candidates.append({'event_index':i,'h_numerator':K.encode(hnum),'h_denominator':K.encode(hden),'h_interval':list(map(str,hbounds)),'h_display_only':list(map(float,hbounds)),'D_numerator':K.encode(numerator),'D_denominator':K.encode(hden),'D_interval':list(map(str,dbounds)),'D_display_only':list(map(float,dbounds)),'strict_positive':dbounds[0]>0,'global_maximizer':i==optimum})
 selected=next(c for c in candidates if c['global_maximizer'])
 out={'status':'POSITIVE_GLOBAL_FREE_H_FIRST_ORDER_DIRECTION_PENDING_REVIEW' if selected['strict_positive'] else ('STRICT_NEGATIVE_GLOBAL_FREE_H_FOR_SELECTED_DISPLACEMENT' if Q(selected['D_interval'][1])<0 else 'UNRESOLVED'),'threshold':'lambda=1+h*epsilon','mass_per_period':'1+epsilon','ratio':'lambda*U/(1+epsilon)','D_formula':'B0+1-C+h*S-sum area*(alpha*h-beta)_+','pure_satellite_halo_kappa':1,'B0_coefficients':K.encode(B),'D_at_h0_coefficients':K.encode(D0),'D_at_h0_display_only':K.disp(D0),'h_event_count':len(events),'weighted_hinge_rows':len(rows),'selected_optimal_event':optimum,'left_h_slope':K.encode(S-cumulative),'right_h_slope':K.encode(S-cumulative-events[optimum][2]),'left_h_slope_display_only':K.disp(S-cumulative),'right_h_slope_display_only':K.disp(S-cumulative-events[optimum][2]),'candidates':candidates,'events':[{'beta':K.encode(b),'alpha':K.encode(a),'weight':K.encode(w),'weighted_intercept':K.encode(t)} for b,a,w,t in events],'normalization_checked':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'free_threshold_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('B0_coefficients','D_at_h0_coefficients','events','candidates','left_h_slope','right_h_slope')},indent=2));print('candidates:',[(c['h_display_only'],c['D_display_only']) for c in candidates])
if __name__=='__main__':
 if sys.argv[1]=='geometry':geometry()
 elif sys.argv[1]=='reflection':run_reflection(int(sys.argv[2]))
 elif sys.argv[1]=='finalize':finalize()
