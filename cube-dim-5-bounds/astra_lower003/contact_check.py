#!/usr/bin/env python3
from pathlib import Path
from functools import lru_cache
from itertools import product
from math import prod
import json,time,hashlib
import threshold_kernel as T
K=T.K;F=K.F;ZERO=K.ZERO;ONE=K.ONE
HERE=Path(__file__).resolve().parent
selected,Z=T.cell();results=json.loads((HERE/'free_threshold_results.json').read_text());candidate=next(x for x in results['candidates'] if x['event_index']==1)
HN=K.decode(candidate['h_numerator']);HD=K.decode(candidate['h_denominator'])
records=[];deficit=ZERO;B0=ZERO
for index in range(18):
 begin=time.monotonic();nodes=0;signs,multiplicity=T.reflection_spec(index);z=tuple(a*s for a,s in zip(Z,signs));saved=json.loads((HERE/f'reflection_{index:02}.json').read_text());B0+=multiplicity*K.decode(saved['outer_intercept']);area_checks=gap_checks=0
 def mask_volume(positive,negative):
  nonlocal_dummy=None
  positive=tuple(sorted(set(positive)));negative=tuple(sorted(set(negative)));axes=[]
  for d in range(4):
   cuts=sorted({x for b in positive+negative for x in b[2*d:2*d+2]});rows=[]
   for lo,hi in zip(cuts,cuts[1:]):
    pm=sum(1<<i for i,b in enumerate(positive) if b[2*d]<=lo and hi<=b[2*d+1])
    nm=sum(1<<i for i,b in enumerate(negative) if b[2*d]<=lo and hi<=b[2*d+1])
    if pm:rows.append((hi-lo,pm,nm))
   axes.append(rows)
  @lru_cache(None)
  def recurse(d,pm,nm):
   global nodes
   nodes+=1
   if nodes>200000 or time.monotonic()-begin>120:raise TimeoutError('registered mask checkpoint')
   if d==4:return ZERO if nm else ONE
   return sum((length*recurse(d+1,pm&p,nm&n) for length,p,n in axes[d] if pm&p),ZERO)
  return recurse(0,(1<<len(positive))-1,(1<<len(negative))-1)
 def count(axis,x,r):return sum(x-r<=integer*K.A[axis]+z[axis] and integer*K.A[axis]+z[axis]<=x+r for integer in range(-4,5))
 for axis,plane,rL,nL,rR,nR,receivers,blockers in T.CONTACTS:
  categories={};other=[j for j in range(5) if j!=axis];faceL=count(axis,plane,rL);faceR=count(axis,plane,rR)
  for b in receivers:
   bands=[]
   for d,j in enumerate(other):
    endpoints={b[2*d],b[2*d+1]}
    for integer in range(-4,5):
     position=integer*K.A[j]+z[j]
     for endpoint in (position-rL,position+rL,position-rR,position+rR):
      if b[2*d]<endpoint and endpoint<b[2*d+1]:endpoints.add(endpoint)
    cuts=sorted(endpoints);bands.append([(lo,hi,count(j,(lo+hi)/2,rL),count(j,(lo+hi)/2,rR)) for lo,hi in zip(cuts,cuts[1:])])
   for pieces in product(*bands):
    key=(faceL*prod(x[2] for x in pieces),faceR*prod(x[3] for x in pieces));categories.setdefault(key,[]).append(tuple(t for x in pieces for t in x[:2]))
  for (sL,sR),boxes in categories.items():
   area=mask_volume(boxes,blockers)
   match=[x for x in saved['hinges'] if x['axis']==axis and x['satellite_counts']==[sL,sR]]
   if area==ZERO:assert not match;continue
   assert len(match)==1 and area==K.decode(match[0]['area']);area_checks+=1
   vl=rL*(HD*sL-HN*nL)/(5*nL);vr=rR*(HD*sR-HN*nR)/(5*nR)
   large=F(10)+(vl if ZERO<vl else -vl)+(vr if ZERO<vr else -vr)
   intervals=((-large,vl),(-vr,large));cuts=sorted({x for b in intervals for x in b})
   union_length=sum((hi-lo for lo,hi in zip(cuts,cuts[1:]) if any(b[0]<=lo and hi<=b[1] for b in intervals)),ZERO)
   loss=2*large-union_length
   alpha=(rL+rR)/5;beta=rL*sL/(5*nL)+rR*sR/(5*nR)
   assert loss==max(ZERO,alpha*HN-beta*HD);gap_checks+=1
   deficit+=multiplicity*area*loss
 records.append({'reflection':index,'nonzero_contact_category_areas_checked':area_checks,'literal_normal_interval_unions_checked':gap_checks,'mask_nodes':nodes,'seconds':time.monotonic()-begin})
 print(json.dumps(records[-1]),flush=True)
geom=json.loads((HERE/'threshold_geometry.json').read_text());C=K.decode(json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients']);KOUT=K.decode(geom['outer_K_full'])
numerator=B0*HD-KOUT*HN-deficit+(HN-HD)*C+HD
assert numerator==K.decode(candidate['D_numerator'])
out={'status':'CONTACT_MASK_AREAS_AND_LITERAL_GAP_UNIONS_VERIFIED','reflections':records,'total_area_checks':sum(r['nonzero_contact_category_areas_checked'] for r in records),'total_gap_checks':sum(r['literal_normal_interval_unions_checked'] for r in records),'normalized_D_numerator_matches':True,'selected_h_event':1,'h_display_only':candidate['h_display_only'],'D_display_only':candidate['D_display_only'],'unfinished_full5D_jet_audit_claimed_passed':False,'shared_geometry_explicit':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'contact_check_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='reflections'},indent=2))
