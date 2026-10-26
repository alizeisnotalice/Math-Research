#!/usr/bin/env python3
"""Read only the saved complete 066-lower and 064-upper fractions."""
from pathlib import Path
from fractions import Fraction
import json,hashlib
p=Path(__file__).resolve().with_name('results.json');d=json.loads(p.read_text())
assert d['status']=='finished' and d['full_union_calls']==0 and not d['new_C6_bound']
lo=Fraction(d['candidate']['exact_interval'][0]);hi=Fraction(d['baseline_complete_064_interval'][1])
assert lo>hi and lo-hi==Fraction(d['candidate_margin_over_complete_064_upper'])
assert len(d['s2_swapped_projection_pairs'])==18
assert all(not z['nonempty'] for z in d['s2_swapped_projection_pairs'])
print(json.dumps({'pass':True,'results_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'new_symbolic_configurations':0,'union_calls':0,'file_writes':0,'new_C6_bound':False}))
