#!/usr/bin/env python3
"""Saved fractions and chain metadata only; no root/field/geometry/union."""
from pathlib import Path
from fractions import Fraction as Q
import json
p=Path(__file__).resolve().parent;r=json.loads((p/'review_results.json').read_text());b=json.loads((p/'frozen_066_baseline.json').read_text());f=r['finite_L1']
assert r['status']=='pass_three_new_hidden_seams_chambers_and_terminal_finite_L1'
old=Q(r['complete_current066_upper']);lo=Q(r['new_ideal_ratio']['exact_interval'][0])
assert old==Q(b['comparison_upper'])==max(Q(v[1]) for v in b['complete_ideal_intervals'].values())
chain=r['new_hidden_chain'];assert len(chain)==3
for i,z in enumerate(chain):
 assert z['all_minmax_choices_unchanged'] and z['next_left_state_equals_current_right']
 assert len(z['projection_pair_certificates'])==[1,16,6][i]
 assert sum(x['required_pair_count'] for x in z['all_tied_function_pairs_and_box_coverage'])==len(z['projection_pair_certificates'])
 assert len(z['new_complete_order_gaps'])==36 and len(z['new_effective_gaps'])==21
 assert all(Q(c['gap']['exact_interval'][0])>0 for c in z['projection_pair_certificates'])
 assert Q(z['new_endpoint_ideal_ratio']['exact_interval'][0])>old
 if i<2:assert z['right_parameter']['coefficients']==chain[i+1]['left_parameter']['coefficients']
assert r['new_hidden_pair_count']==23
assert r['new_ideal_ratio']==chain[-1]['new_endpoint_ideal_ratio']
assert r['final_boundary_ties']==[[31,7],[3,29]] and not r['final_right_state_or_slope_classification_performed'] and r['final_projection_checks']==0
eps=Q(f['epsilon']);theta=Q(f['threshold']);K=f['receiver_index_radius']
assert eps==Q(1,10**22) and Q(f['bump_halfwidth'])==eps/2
assert theta==1-Q(1,10**16) and K==10**16 and f['source_index_radius']==K+2
loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor;target=Q(2071938494109,10**12)
assert loss==Q(f['union_loss_upper']) and inside==Q(f['inward_ratio_lower']) and factor==Q(f['boundary_factor'])
assert finite==Q(f['weak_ratio_lower']) and target==Q(f['strict_rational_target']) and finite>target>old
assert finite-target==Q(f['margin_over_target']) and finite-old==Q(f['margin_over_complete066_upper'])
assert lo-old==Q(r['new_ideal_margin_over_complete066_upper'])
print('PASS saved exact hidden-chain metadata, finite L1, and complete066-upper comparisons; zero new scientific configurations')
