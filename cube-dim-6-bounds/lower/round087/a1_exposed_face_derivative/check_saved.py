#!/usr/bin/env python3
"""Saved structural/cap consistency only; no DP, field, or union execution."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,resource,signal
resource.setrlimit(resource.RLIMIT_CPU,(20,21));signal.alarm(180)
h=Path(__file__).resolve().parent;a=json.loads((h/'results.json').read_text());d=json.loads((h/'structure_results.json').read_text())
assert a['status']=='stopped_exception' and a['error']=="RuntimeError('DP/configuration cap')"
assert a['new_DP_nodes']==1001 and a['sides']=={}
assert d['prior_DP_admitted_nodes']==1000 and d['prior_DP_counter_including_rejected_request']==1001
assert d['status']=='finished_structure_only' and d['new_DP_nodes']==0 and d['new_5D_union_calls']==0
assert a['new_configurations']+d['new_configurations']==3
assert len(a['normal_functions'])==24 and len(a['tail_rank_coordinates'])==21
assert d['minmax_selections_unchanged_across_current']
for name,tie,count in [('left',[[13,6]],13),('right',[[6,13]],15)]:
 s=d['sides'][name];assert s['positive_boxes']==127 and s['all_tie_groups']==tie
 assert len(s['required_weighted_slabs'])==count
 assert all(x['5D_union_area']=='NOT_COMPUTED' for x in s['required_weighted_slabs'])
 assert Fraction(s['trivial_derivative_lower']['exact_interval'][0])<0<Fraction(s['trivial_derivative_upper']['exact_interval'][1])
assert d['input_sha256'][str(h/'results.json')]==hashlib.sha256((h/'results.json').read_bytes()).hexdigest()
print(json.dumps({'pass':True,'structure_results_sha256':hashlib.sha256((h/'structure_results.json').read_bytes()).hexdigest(),'new_configurations':0,'new_DP_nodes':0,'union_calls':0,'writes':0,'derivative_sign':'UNRESOLVED'}));signal.alarm(0)
