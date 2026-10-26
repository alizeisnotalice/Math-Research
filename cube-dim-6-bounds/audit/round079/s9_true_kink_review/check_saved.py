#!/usr/bin/env python3
"""Portable proof-table checker: no author imports, new search, or union solver."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,time
T=time.process_time();W=time.monotonic()
resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def timeout(*_):raise TimeoutError('180-second proof-table check')
signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
p=Path(__file__).resolve().parent/'review_results.json'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='f919dd766ca72851bf48099c55c3b38df0e012b4c8a3d334e1312f6bbc5d9615'
r=json.loads(p.read_text());assert r['status']=='pass_complete_s9_kink_and_local_maximum'
def parse(d):return {tuple(map(int,k.split(','))):Q(v) for k,v in d.items() if Q(v)}
def add(a,b):
 c=dict(a)
 for k,v in b.items():c[k]=c.get(k,Q(0))+v
 return {k:v for k,v in c.items() if v}
def scale(a,q):return {k:v*q for k,v in a.items() if v*q}
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
 c={}
 for (i,j),v in a.items():
  for (k,l),w in b.items():
   ij=(i+k)%6,(j+l)%6;c[ij]=c.get(ij,Q(0))+v*w*2**((i+k)//6)*3**((j+l)//6)
 return {k:v for k,v in c.items() if v}
def const(q):return {(0,0):Q(q)} if q else {}
DEN=10**100
def root(n):
 lo,hi=DEN,2*DEN;target=n*DEN**6
 while hi-lo>1:
  m=(lo+hi)//2
  if m**6<=target:lo=m
  else:hi=m
 assert lo**6<=target<(lo+1)**6
 return Q(lo,DEN),Q(lo+1,DEN)
ul,uh=root(2);vl,vh=root(3)
powers={(i,j):(ul**i*vl**j,uh**i*vh**j) for i in range(6) for j in range(6)}
def bounds(a):
 lo=hi=Q(0)
 for k,v in a.items():
  l,h=powers[k];lo+=v*(l if v>0 else h);hi+=v*(h if v>0 else l)
 return lo,hi
def compare(a,b):
 d=sub(a,b)
 if not d:return 0
 lo,hi=bounds(d);assert lo>0 or hi<0
 return 1 if lo>0 else -1
def view(v):
 a=parse(v['coefficients']);lo,hi=bounds(a)
 assert [str(lo),str(hi)]==v['exact_interval'];return a
for k in ['s8','s9','s10','exact_background_5D_union','exact_uncovered_5D_volume','left_slope','exact_ratio_slope_jump','right_slope']:view(r[k])
cert=r['union_certificate'];coords=[view(x) for x in cert['coordinate_catalog']]
assert all(compare(a,b)<0 for a,b in zip(coords,coords[1:]))
S=cert['states'];assert len(S)==56 and [x['id'] for x in S]==list(range(len(S)))
boxkeys=[tuple(tuple(tuple(p) for p in B) for B in x['boxes']) for x in S]
for i,st in enumerate(S):
 boxes=boxkeys[i];dim=st['dimension'];assert boxes==tuple(sorted(set(boxes))) and boxes
 assert all(len(B)==dim and all(0<=l<h<len(coords) for l,h in B) for B in boxes)
 if dim==0:
  assert boxes==((),) and st['edges']==[] and parse(st['volume'])==const(1);continue
 cuts=sorted({q for B in boxes for q in B[0]});expected=[]
 for l,h in zip(cuts,cuts[1:]):
  child=tuple(sorted({B[1:] for B in boxes if B[0][0]<=l and h<=B[0][1]}))
  if child:expected.append((l,h,child))
 assert len(expected)==len(st['edges']);value={}
 for (l,h,child),edge in zip(expected,st['edges']):
  assert edge['lo']==l and edge['hi']==h
  j=edge['child'];assert S[j]['dimension']==dim-1 and boxkeys[j]==child
  value=add(value,mul(sub(coords[h],coords[l]),parse(S[j]['volume'])))
 assert value==parse(st['volume'])
# Root corresponds exactly to the saved 21 genuine background sections.
def coordinate_id(a):return coords.index(parse(a))
retained=r['retained_cross_sections'];assert len(retained)==21
root=tuple(sorted(tuple((coordinate_id(a),coordinate_id(b)) for a,b in reversed(B)) for B in retained))
assert root==boxkeys[cert['root_state']] and cert['original_axis_order']==[6,5,4,3,1]
covered=parse(S[cert['root_state']]['volume']);assert covered==parse(r['exact_background_5D_union']['coefficients'])
uncovered=sub(const(Q(1,32)),covered);jump=scale(uncovered,-64);left=parse(r['left_slope']['coefficients']);right=add(left,jump)
assert uncovered==parse(r['exact_uncovered_5D_volume']['coefficients']) and jump==parse(r['exact_ratio_slope_jump']['coefficients']) and right==parse(r['right_slope']['coefficients'])
assert compare(uncovered,{})>0 and compare(left,{})>0 and compare(right,{})<0
records={x['raw_index']:x for x in r['all125_background_records']};assert set(records)==set(range(127))-{0,16}
assert {c:sum(x['classification']==c for x in records.values()) for c in ('active','below','above')}==r['background_counts']=={'active':76,'below':43,'above':6}
positive={k:x for k,x in records.items() if x.get('positive_5D')};assert len(positive)==44
proof=r['all44_contained_in_retained21'];assert {x['raw_index'] for x in proof}==set(positive) and len(proof)==44
for x in proof:
 B=positive[x['raw_index']]['cross_section'];C=retained[x['retained_index']]
 assert all(compare(parse(c),parse(a))<=0 and compare(parse(b),parse(d))<=0 for (a,b),(c,d) in zip(B,C))
assert all(B in [x['cross_section'] for x in positive.values()] for B in retained)
for x in records.values():
 for v in x['normal_margins']:assert compare(view(v),{})>0
for x in r['normal_collar_endpoint_certificates']:assert x['all125_pass'] and compare(view(x['minimum_strict_background_strip_margin']),{})>0
assert Q(r['normal_collar_radius'])>0 and len(r['normal_collar_endpoint_certificates'])==2
assert len(r['new_right_chamber_gap_certificate'])==36
for x in r['new_right_chamber_gap_certificate']:
 a=parse(x['left_gap']);b=parse(x['right_gap']);assert compare(a,{})>=0 and compare(b,{})>=0 and (a or b)
assert r['unique_simultaneous_tie']==[7,2] and r['next_ties']==[[33,28]]
assert r['origin_covered_by_retained_indices'] and r['open_uncovered_box']['not_used_as_full_jump']
assert not r['new_C6_bound'] and not r['global_optimum_claim']
assert r['resource']['configurations']==7 and r['resource']['local_nodes']==786 and r['resource']['new_5D_union_states']==56 and r['resource']['full_6D_union_calls']==0
print(json.dumps({'status':'PASS','certificate_states_verified':56,'new_union_searches':0,'new_scientific_configurations':0,'source_result_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'CPU_seconds':time.process_time()-T,'wall_seconds':time.monotonic()-W},indent=2))
signal.alarm(0)
