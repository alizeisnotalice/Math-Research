#!/usr/bin/env python3
"""Saved-only Fraction/vector consistency; no field or union execution."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,resource,signal
resource.setrlimit(resource.RLIMIT_CPU,(20,21));signal.alarm(180)
h=Path(__file__).resolve().parent;s=json.loads((h/'results.json').read_text());d=json.loads((h/'section_results.json').read_text())
assert s['status']=='stopped_exception' and s['error']=='AssertionError()' and s['new_5D_union_nodes']==0
assert d['status']=='finished_strict_local_maximum' and d['full_6D_union_calls']==0 and not d['new_C6_bound']
assert d['new_5D_union_calls']==1 and d['old_phase_reruns']==0 and d['old_background_nodes_replayed']==0
assert d['combined_nodes']<=1000 and d['combined_symbolic_configurations']<=8
assert len(s['all_background_normal_classifications'])==125
for z in s['all_background_normal_classifications']:
 if z['classification']=='strictly_active_at_height':
  assert Q(z['strict_normal_collar_below']['exact_interval'][0])>0 and Q(z['strict_normal_collar_above']['exact_interval'][0])>0
 else:assert Q(z['strict_normal_distance']['exact_interval'][0])>0
assert Q(s['explicit_rational_common_normal_neighborhood'])>0
assert Q(s['left_slope']['exact_interval'][0])>0 and Q(d['right_slope']['exact_interval'][1])<0
assert Q(d['exact_uncovered_5D_volume']['exact_interval'][0])>0
keys=[f'{i},{j}' for i in range(6) for j in range(6)]
def c(x,k):return Q(x['coefficients'].get(k,'0'))
for k in keys:
 p=Q(1,32) if k=='0,0' else Q(0)
 assert c(d['exact_background_5D_union'],k)+c(d['exact_uncovered_5D_volume'],k)==p
 assert c(d['exact_ratio_slope_jump'],k)==-64*c(d['exact_uncovered_5D_volume'],k)
 assert c(d['right_slope'],k)==c(s['left_slope'],k)+c(d['exact_ratio_slope_jump'],k)
for a,b in d['open_uncovered_region']['rational_inward_open_intervals']:assert Q(a)<Q(b)
assert len(d['open_uncovered_region']['all_21_retained_boxes_disjoint_witnesses'])==21
assert d['next_states']['left']==s['current_right_state_reused']
assert d['input_sha256'][str(h/'results.json')]==hashlib.sha256((h/'results.json').read_bytes()).hexdigest()
print(json.dumps({'pass':True,'section_results_sha256':hashlib.sha256((h/'section_results.json').read_bytes()).hexdigest(),'new_symbolic_configurations':0,'union_calls':0,'writes':0,'new_C6_bound':False}))
signal.alarm(0)
