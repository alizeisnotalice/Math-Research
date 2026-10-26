#!/usr/bin/env python3
"""Saved rational source-count/root and conservative-volume verification."""
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import json,hashlib,resource,signal
resource.setrlimit(resource.RLIMIT_CPU,(20,21));signal.alarm(180)
p=Path(__file__).resolve().with_name('results.json');d=json.loads(p.read_text())
assert d['status']=='finished' and d['union_calls']==0 and d['old_program_replays']==0 and not d['new_C6_bound']
oldL,oldU=map(Q,d['old_full_075_ideal_interval']);cL,cU=map(Q,d['period_cell_volume_interval'])
assert cU/2<oldL and Q(d['certified_positive_headroom'])==cL/oldU-1
w=Q(d['derived_weight']);headroom=Q(d['certified_positive_headroom']);assert w==Q(1,16) and w<=headroom/2<2*w
assert len(d['candidate_witnesses'])==3

def vol(B):
 if any(b<=a for a,b in B):return Q(0)
 return prod(b-a for a,b in B)
def box(z):return [(Q(a),Q(b)) for a,b in z]
total=Q(0)
for z in d['candidate_witnesses']:
 m=prod(b-a+1 for a,b in zip(z['unit_lo'],z['unit_hi']));k=prod(b-a+1 for a,b in zip(z['sat_lo'],z['sat_hi']))
 assert m==z['unit_count'] and k==z['satellite_count']
 M=Q(z['selected_mass']);assert M==m+k*w
 lo,hi=Q(z['radius_lower']),Q(z['radius_upper']);assert lo**6<=M/64<hi**6
 for i,a in enumerate(z['true_source_extrema']):
  values=[Q(z['unit_lo'][i]),Q(z['unit_hi'][i]),Q(z['sat_lo'][i])+Q(1,2),Q(z['sat_hi'][i])+Q(1,2)]
  assert [str(min(values)),str(max(values))]==a['actual_combined_normalized_coordinate_extrema']
 inner,outer,ref=map(box,[z['rational_inward_box'],z['rational_outer_box'],z['reference_old_inner_box']])
 assert all(a<c<e<b or a<=c<e<=b for (a,b),(c,e) in zip(outer,inner))
 intersection=[(max(a,c),min(b,e)) for (a,b),(c,e) in zip(outer,ref)]
 du=vol(outer)-vol(intersection);assert du==Q(z['outer_minus_reference_inner_volume_upper']);total+=du
assert total==Q(d['sum_new_volume_increment_upper'])
upper=(oldU+64*total)/(1+w);assert upper==Q(d['full_three_mixed_box_family_ratio_upper']) and upper<oldL
print(json.dumps({'pass':True,'results_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'new_geometric_configurations':0,'union_calls':0,'writes':0,'new_C6_bound':False}))
signal.alarm(0)
