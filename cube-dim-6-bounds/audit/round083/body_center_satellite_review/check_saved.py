#!/usr/bin/env python3
"""Portable exact verification of the saved E083 finite certificates; no union."""
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import hashlib,json,resource,signal,time
T=time.process_time();W=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def timeout(*_):raise TimeoutError('180 seconds including cleanup')
signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
p=Path(__file__).resolve().with_name('review_results.json')
assert hashlib.sha256(p.read_bytes()).hexdigest()=='6099208470aa234d10f4f34e6665e8f4181e78fbd9ee94c739ebae2ec2b550db'
r=json.loads(p.read_text());assert r['status']=='pass_explicit_three_mixed_subfamily_upper'
DEN=10**100
def root(n):
 lo,hi=DEN,2*DEN
 while hi-lo>1:
  mid=(lo+hi)//2
  if mid**6<=n*DEN**6:lo=mid
  else:hi=mid
 assert lo**6<=n*DEN**6<(lo+1)**6
 return Q(lo,DEN),Q(lo+1,DEN)
ul,uh=root(2);vl,vh=root(3)
def field_bounds(c):
 lo=hi=Q(0)
 for ij,value in c.items():
  i,j=map(int,ij.split(','));q=Q(value);l,h=ul**i*vl**j,uh**i*vh**j
  lo+=q*(l if q>=0 else h);hi+=q*(h if q>=0 else l)
 return lo,hi
def view(v):
 b=field_bounds(v['coefficients']);assert b==tuple(map(Q,v['exact_interval']));return b
def box(b):return [tuple(map(Q,p)) for p in b]
def volume(b):return prod(hi-lo for lo,hi in b) if all(lo<hi for lo,hi in b) else Q(0)
ab=[view(x) for x in r['spacing']];cell=view(r['period_cell_volume']);dcL,dcU=map(Q,r['author_period_cell_interval']);assert dcL<=cell[0]<=cell[1]<=dcU
assert r['spacing'][0]['coefficients']=={'0,0':'1','4,0':'1/2'} and all(x['coefficients']=={'0,0':'1/2','1,0':'1/2'} for x in r['spacing'][1:])
oldL,oldU=Q(r['comparison_lower']),Q(r['comparison_upper']);assert oldL==min(Q(x[0]) for x in r['complete_E075_intervals'].values()) and oldU==max(Q(x[1]) for x in r['complete_E075_intervals'].values())
w=Q(r['weight']);assert w==Q(1,16) and Q(r['cell_mass'])==1+w
headroom=dcL/oldU-1;assert headroom==Q(r['positive_headroom']) and w<=headroom/2<2*w
assert Q(r['equal_weight_level1_ceiling_upper'])==dcU/2<oldL
total=Q(0);axis_count=0
for z in r['witnesses']:
 Ulo,Uhi=z['unit_index_box'];Slo,Shi=z['satellite_integer_index_box']
 m=prod(b-a+1 for a,b in zip(Ulo,Uhi));k=prod(b-a+1 for a,b in zip(Slo,Shi));M=Q(z['mass'])
 assert m==z['unit_count'] and k==z['satellite_count'] and M==m+k*w
 lo,hi=map(Q,z['critical_radius_enclosure']);eta=Q(z['strict_source_margin']);rt=Q(z['one_common_strict_test_radius'])
 assert 0<lo<hi and lo**6<=M/64<hi**6 and rt==lo-eta/2>0 and (2*rt)**6<M
 assert list(map(Q,z['root_sixth_power_residuals']))==[M/64-lo**6,hi**6-M/64] and Q(z['strict_level1_mass_residual'])==M-(2*rt)**6
 inn,outer,ref=map(box,[z['rational_inward_box'],z['rational_outer_box'],z['reference_old_inner_box']])
 for j,x in enumerate(z['axis_certificates']):
  axis_count+=1;assert x['axis']==j+1
  endpoints=[Q(Ulo[j]),Q(Uhi[j]),Q(Slo[j])+Q(1,2),Q(Shi[j])+Q(1,2)];mn,mx=min(endpoints),max(endpoints)
  assert list(map(Q,x['normalized_source_extrema']))==[mn,mx]
  def coord(q):return (q*ab[j][0],q*ab[j][1]) if q>=0 else (q*ab[j][1],q*ab[j][0])
  mnB,mxB=coord(mn),coord(mx);assert mnB==tuple(map(Q,x['physical_source_min_interval'])) and mxB==tuple(map(Q,x['physical_source_max_interval']))
  il,ih=inn[j];ol,oh=outer[j];rl,rh=ref[j]
  assert 0<=il<ih<=ab[j][0]/2 and ol<=il<ih<=oh
  slacks=[il-(mxB[1]-lo+eta),(mnB[0]+lo-eta)-ih]
  assert slacks==list(map(Q,x['aggregate_source_slacks'])) and min(slacks)>=0
  assert ol<=max(Q(0),mxB[0]-hi) and oh>=min(ab[j][1]/2,mnB[1]+hi)
  oldB=[field_bounds(c) for c in x['old_reference_exact_endpoints']]
  assert oldB==[tuple(map(Q,b)) for b in x['old_reference_endpoint_intervals']]
  assert oldB[0][1]+eta<=rl<rh<=oldB[1][0]-eta
  assert x['unit_reflected_range']==[-Uhi[j],-Ulo[j]] and x['satellite_reflected_integer_range']==[-Shi[j]-1,-Slo[j]-1]
 inter=[(max(a,c),min(b,d)) for (a,b),(c,d) in zip(outer,ref)];assert inter==box(z['outer_intersection_reference_inner'])
 du=volume(outer)-volume(inter);assert du==Q(z['increment_volume_upper'])>=0 and volume(outer)==Q(z['outer_volume']) and volume(inter)==Q(z['intersection_volume']) and volume(inn)>0
 total+=du
assert axis_count==18 and len(r['witnesses'])==3 and total==Q(r['sum_increment_upper'])
upper=(oldU+64*total)/(1+w);assert upper==Q(r['family_ratio_upper'])<oldL and oldL-upper==Q(r['margin_below_complete_E075_lower'])
note=r['normalization_threshold_note'];assert Q(note['necessary_threshold_lower'])==w*oldL/64 and Q(note['threshold_upper'])==w*oldU/64 and Q(note['author_saved_safe_sufficient_threshold'])==((1+w)*oldU-oldL)/64
assert not r['actual_union_volume_computed'] and not r['new_C6_bound'] and r['resource']['old144_source_endpoint_nodes_replayed']==r['resource']['union_calls']==0
print(json.dumps({'status':'PASS','saved_axis_certificates_verified':18,'root_endpoint_sixth_power_pairs_verified':3,'old144_endpoint_nodes_replayed':0,'union_calls':0,'new_geometric_configurations':0,'CPU_seconds':time.process_time()-T,'wall_seconds':time.monotonic()-W},indent=2))
signal.alarm(0)
