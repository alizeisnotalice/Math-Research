#!/usr/bin/env python3
"""Saved-only exact endpoint/proof-chain sanity; no field or union imports."""
from pathlib import Path
from fractions import Fraction
import json,hashlib
p=Path(__file__).resolve().with_name('results.json');d=json.loads(p.read_text())
assert d['status']=='stopped_overlap' and d['full_union_calls']==0 and not d['new_C6_bound']
assert d['event_classifications']==4 and len(d['events'])==4
assert d['symbolic_configurations']<=32 and d['local_nodes']<=4000
old=Fraction(d['baseline_complete_066_interval'][1])
for i,e in enumerate(d['events'][:-1]):
 assert e['classification']=='hidden_strict_fixed_projection_disjointness'
 assert e['left_state']['selections']==e['right_state']['selections']
 assert all(z['strictly_disjoint'] and Fraction(z['gap']['exact_interval'][0])>0 for z in e['projection_pair_checks'])
 nxt=e['next_interval'];assert nxt['right_parameter']['coefficients']==d['events'][i+1]['parameter']['coefficients']
 c=nxt['candidate'];assert Fraction(c['ratio']['exact_interval'][0])>old
 assert Fraction(c['ratio']['exact_interval'][0])-old==Fraction(c['exact_margin_over_complete_066_upper'])
assert d['best_new_candidate']==d['events'][-2]['next_interval']['candidate']
stop=d['events'][-1];assert stop['classification']=='stop_first_nonempty_or_touching_projection'
assert stop['projection_pair_checks'][-1]['all_five_intervals_positive'] is True
print(json.dumps({'pass':True,'results_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'new_symbolic_configurations':0,'union_calls':0,'writes':0,'new_C6_bound':False}))
