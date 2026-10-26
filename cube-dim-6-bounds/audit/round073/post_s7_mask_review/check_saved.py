#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as Q
import json,resource,signal
resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180 seconds including cleanup')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
p=Path(__file__).resolve().parent;r=json.loads((p/'review_results.json').read_text());b=json.loads((p/'frozen_070_baseline.json').read_text());f=r['finite_L1']
assert r['status']=='pass_new_s7_disjoint_seam_chamber_and_s8_finite_L1'
assert r['all_simultaneous_ties']==[[11,26]] and r['new_s7_right_germ_verified'] and r['all_127_minmax_choices_unchanged']
assert len(r['all_seven_tie_labels_effective'])==7 and r['complete_effective_groups']==[[31],[81,84,86,93,95,100]]
assert len(r['six_disjoint_projection_pairs'])==6 and all(Q(x['gap']['exact_interval'][0])>0 for x in r['six_disjoint_projection_pairs'])
assert r['existing_cover_candidates_tested']==0 and len(r['new_complete_order_gaps'])==36 and len(r['new_effective_gaps'])==21
assert r['next_ties']==[[31,2]] and r['s8_left_state_equals_s7_right']
old=Q(r['complete_current070_upper']);lo=Q(r['new_ideal_ratio']['exact_interval'][0]);assert old==Q(b['comparison_upper'])==max(Q(v[1]) for v in b['complete_ideal_intervals'].values())
eps=Q(f['epsilon']);theta=Q(f['threshold']);K=f['receiver_index_radius'];assert eps==Q(1,10**22) and Q(f['bump_halfwidth'])==eps/2 and theta==1-Q(1,10**16) and K==10**16 and f['source_index_radius']==K+2
loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor;target=Q(2072477059470,10**12)
assert loss==Q(f['union_loss_upper']) and inside==Q(f['inward_ratio_lower']) and factor==Q(f['boundary_factor'])
assert finite==Q(f['weak_ratio_lower']) and target==Q(f['strict_rational_target']) and finite>target>old
assert finite-target==Q(f['margin_over_target']) and finite-old==Q(f['margin_over_complete070_upper']) and lo-old==Q(r['new_ideal_margin_over_complete070_upper'])
print('PASS saved six-pair completeness, finite L1 and full070 comparisons; zero new scientific configurations',flush=True)
signal.alarm(0)
