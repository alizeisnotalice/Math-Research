#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as Q
import json
p=Path(__file__).resolve().parent;r=json.loads((p/'review_results.json').read_text());b=json.loads((p/'frozen_068_baseline.json').read_text());f=r['finite_L1']
assert r['status']=='pass_complete_cover_mask_new_interval_and_finite_L1'
assert r['all_127_minmax_choices_unchanged'] and r['new_s6_right_germ_verified']
assert len(r['six_complete_cover_pairs'])==6 and len(r['eight_inactive_labels'])==8
for c in r['six_complete_cover_pairs']:
 assert c['cover_mass']==8
 for face in c['fixed_inclusion_margins']:
  assert Q(face['low_gap']['exact_interval'][0])>=0 and Q(face['high_gap']['exact_interval'][0])>=0
 assert Q(c['normal_collar_below']['exact_interval'][0])>0 and Q(c['normal_collar_above']['exact_interval'][0])>0
collar=r['two_sided_strip_collar'];assert Q(collar['parameter_radius'])==Q(1,10**6) and collar['inside_preceding_and_new_certified_chambers']
assert [x['side'] for x in collar['endpoint_checks']]==[-1,1]
for e in collar['endpoint_checks']:
 assert len(e['cover_checks'])==6 and len(e['inactive_checks'])==8
 for c in e['cover_checks']:
  for z in c['both_swapped_endpoints_strictly_inside']:
   assert Q(z['from_lower']['exact_interval'][0])>0 and Q(z['from_upper']['exact_interval'][0])>0
 assert all(Q(z['gap']['exact_interval'][0])>0 for z in e['inactive_checks'])
old=Q(r['complete_current068_upper']);lo=Q(r['new_ideal_ratio']['exact_interval'][0]);assert old==Q(b['comparison_upper'])==max(Q(v[1]) for v in b['complete_ideal_intervals'].values())
eps=Q(f['epsilon']);theta=Q(f['threshold']);K=f['receiver_index_radius'];assert eps==Q(1,10**22) and Q(f['bump_halfwidth'])==eps/2 and theta==1-Q(1,10**16) and K==10**16 and f['source_index_radius']==K+2
loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor;target=Q(2071993199356,10**12)
assert loss==Q(f['union_loss_upper']) and inside==Q(f['inward_ratio_lower']) and factor==Q(f['boundary_factor'])
assert finite==Q(f['weak_ratio_lower']) and target==Q(f['strict_rational_target']) and finite>target>old
assert finite-target==Q(f['margin_over_target']) and finite-old==Q(f['margin_over_complete068_upper']) and lo-old==Q(r['new_ideal_margin_over_complete068_upper'])
print('PASS saved complete two-sided strip margins, finite L1 and full068 comparisons; zero new scientific configurations')
