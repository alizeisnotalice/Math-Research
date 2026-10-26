#!/usr/bin/env python3
"""Read-only exact comparison of the two saved full ideal endpoints."""
from pathlib import Path
from fractions import Fraction
import json,hashlib
path=Path(__file__).resolve().with_name('results.json')
d=json.loads(path.read_text())
assert d['status']=='finished' and d['full_union_calls']==0 and not d['new_C6_bound']
new_lower=Fraction(d['candidate']['exact_interval'][0])
old_upper=Fraction(d['baseline_058_ideal_interval'][1])
assert new_lower>old_upper
assert new_lower-old_upper==Fraction(d['candidate_margin_over_complete_058_upper'])
print(json.dumps({'pass':True,'saved_result_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'new_configurations':0,'union_calls':0,'new_C6_bound':False}))
