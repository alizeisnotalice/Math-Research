#!/usr/bin/env python3
"""Saved-only rational certificate sanity; no geometry/field/union import."""
from pathlib import Path
from fractions import Fraction
import json,hashlib
p=Path(__file__).resolve().with_name('results.json');d=json.loads(p.read_text())
assert d['status']=='finished' and d['full_union_calls']==0 and not d['new_C6_bound']
assert d['all_simultaneous_effective_strips_paid'] and d['seam_hidden_and_slope_preserved']
first,second=d['tie_diagnostics'];assert first['ids']==[31,7] and second['ids']==[3,29]
assert len(first['pair_covers'])==6 and len(second['checks'])==8
for c in first['pair_covers']:
 assert c['covered_with_strict_relative_cell_collar'] and len(c['fixed_face_checks'])==5
 assert Fraction(c['normal_collar_below']['exact_interval'][0])>0
 assert Fraction(c['normal_collar_above']['exact_interval'][0])>0
 for q in c['fixed_face_checks']:
  assert Fraction(q['low_gap']['exact_interval'][0])>0 or q['low_equality_is_cell_face']
  assert Fraction(q['high_gap']['exact_interval'][0])>0 or q['high_equality_is_cell_face']
for c in second['checks']:assert Fraction(c['unused_upper_minus_halfcell']['exact_interval'][0])>0
lo=Fraction(d['candidate']['exact_interval'][0]);hi=Fraction(d['baseline_complete_E068_interval'][1])
assert lo>hi and lo-hi==Fraction(d['candidate_margin_over_complete_E068_upper'])
assert d['next_states']['left']==d['current_right_state_reused']
print(json.dumps({'pass':True,'results_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'new_symbolic_configurations':0,'union_calls':0,'writes':0,'new_C6_bound':False}))
