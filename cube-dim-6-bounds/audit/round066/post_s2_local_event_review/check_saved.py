#!/usr/bin/env python3
"""Exact saved-Fraction checks only: no roots, field, geometry or union."""
from pathlib import Path
from fractions import Fraction as Q
import json
p=Path(__file__).resolve().parent
r=json.loads((p/'review_results.json').read_text());b=json.loads((p/'frozen_064_baseline.json').read_text());f=r['finite_L1']
assert r['status']=='pass_new_s2_seam_new_interval_and_finite_L1'
lo=Q(r['new_ideal_ratio']['exact_interval'][0]);old=Q(r['complete_current064_upper'])
assert old==Q(b['comparison_upper'])==max(Q(v[1]) for v in b['complete_ideal_intervals'].values())
eps=Q(f['epsilon']);theta=Q(f['threshold']);K=f['receiver_index_radius']
assert eps==Q(1,10**22) and Q(f['bump_halfwidth'])==eps/2
assert theta==1-Q(1,10**16) and K==10**16 and f['source_index_radius']==K+2
loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor;target=Q(2070186602851,10**12)
assert loss==Q(f['union_loss_upper']) and inside==Q(f['inward_ratio_lower'])
assert factor==Q(f['boundary_factor']) and finite==Q(f['weak_ratio_lower'])
assert target==Q(f['strict_rational_target']) and finite>target>old
assert finite-target==Q(f['margin_over_target']) and finite-old==Q(f['margin_over_complete064_upper'])
assert lo-old==Q(r['new_ideal_margin_over_complete064_upper'])
assert len(r['s2_projection_pair_certificate'])==18 and r['s2_hidden_by_complete_crossing_pair_disjointness']
print('PASS saved exact finite-L1 and complete064-upper comparisons; zero new geometry/field/union')
