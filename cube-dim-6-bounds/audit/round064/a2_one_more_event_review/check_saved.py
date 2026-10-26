#!/usr/bin/env python3
"""Saved rational certificate check only; no geometry, field or union."""
from pathlib import Path
from fractions import Fraction as Q
import json
p=Path(__file__).resolve().parent
r=json.loads((p/'review_results.json').read_text());f=r['finite_L1']
assert r['status']=='pass_new_interval_and_finite_L1'
lo=Q(r['new_ideal_ratio']['exact_interval'][0]);old=Q(r['complete_current058_upper'])
s=json.loads((p/'frozen_lower_ledger.json').read_text())
assert old==Q(s['lower']['exact_interval'][1])
eps=Q(f['epsilon']);theta=Q(f['threshold']);K=f['receiver_index_radius']
assert eps==Q(1,10**22) and Q(f['bump_halfwidth'])==eps/2
assert theta==1-Q(1,10**16) and K==10**16 and f['source_index_radius']==K+2
loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6
finite=theta*inside*factor;target=Q(2068316441022,10**12)
assert loss==Q(f['union_loss_upper']) and inside==Q(f['inward_ratio_lower'])
assert factor==Q(f['boundary_factor']) and finite==Q(f['weak_ratio_lower'])
assert target==Q(f['strict_rational_target']) and finite>target>old
assert finite-target==Q(f['margin_over_target'])
assert finite-old==Q(f['margin_over_complete058_upper'])
assert lo-old==Q(r['new_ideal_margin_over_complete058_upper'])
print('PASS saved-only exact finite-L1 and complete058-upper comparisons; zero new geometry/field/union')
