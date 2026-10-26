#!/usr/bin/env python3
"""Saved-only rational comparison. No field import, geometry, or union."""
from fractions import Fraction
from pathlib import Path
import json
import hashlib

HERE=Path(__file__).resolve().parent
path=HERE/'results.json'
d=json.loads(path.read_text())
candidate_lower=Fraction(d['next_affine_ideal_candidate']['exact_interval'][0])
old008_upper=Fraction(d['current_candidate_interval'][1])
margin=candidate_lower-old008_upper
assert margin>0
assert margin==Fraction(d['candidate_margin_over_complete_current_upper'])
assert d['exact_above_complete_current_ideal_upper'] is True
assert d['full_union_calls']==0 and d['new_C6_bound'] is False
print(json.dumps({
 'saved_results_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
 'candidate_lower_field':'next_affine_ideal_candidate.exact_interval[0]',
 'old008_upper_field':'current_candidate_interval[1]',
 'exact_fraction_comparison_passed':True,
 'stored_margin_exactly_matched':True,
 'new_symbolic_configurations':0,'union_calls':0,'new_C6_bound':False
},indent=2))
