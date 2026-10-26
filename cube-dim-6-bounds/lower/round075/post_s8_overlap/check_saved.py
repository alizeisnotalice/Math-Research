#!/usr/bin/env python3
"""Saved-only rational checks; no geometry, field, or union execution."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,signal,resource
resource.setrlimit(resource.RLIMIT_CPU,(20,21));signal.alarm(180)
p=Path(__file__).resolve().with_name('results.json');d=json.loads(p.read_text())
assert d['status']=='finished' and d['full_union_calls']==0 and not d['new_C6_bound']
assert d['tie_edges']==[[31,2]] and len(d['tie_diagnostics'])==1
r=d['tie_diagnostics'][0]
assert len(r['all_labels'])==7 and all(z['function_is_effectively_selected'] for z in r['all_labels'])
assert len(r['pair_diagnostics'])==6
assert all(z['classification']=='strictly_disjoint' and Fraction(z['positive_gap']['exact_interval'][0])>0 for z in r['pair_diagnostics'])
assert d['existing_cover_candidates_tested']==0 and d['seam_hidden_and_g_preserved']
assert d['next_states']['left']==d['current_right_state_reused']
lo=Fraction(d['candidate']['exact_interval'][0]);hi=Fraction(d['baseline_complete_E073_interval'][1])
assert lo>hi and lo-hi==Fraction(d['candidate_margin_over_complete_E073_upper'])
print(json.dumps({'pass':True,'results_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'new_symbolic_configurations':0,'union_calls':0,'writes':0,'new_C6_bound':False}))
signal.alarm(0)
