#!/usr/bin/env python3
"""Exact dynamic-ledger comparison for the two inward rational certificates."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/audit/lower003')
HERE = Path(__file__).resolve().parent
k1 = json.loads((MAIN / 'rational_inward_candidate_results.json').read_text())
k16 = json.loads((MAIN / 'rational_inward_k16_results.json').read_text())
r1 = Fraction(k1['exact_rational_lower'])
r16 = Fraction(k16['exact_rational_lower'])
t16 = Fraction(k16['simple_strict_lower'])
result = {
    'status': 'pass',
    'comparison': 'current k16 exact rational inward union versus prior k1 exact rational inward union',
    'k16_exact_rational_lower_strictly_exceeds_k1_exact_rational_lower': r16 > r1,
    'k16_finite_L1_target_strictly_exceeds_k1_exact_rational_lower': t16 > r1,
    'shared_denominator': r16.denominator == r1.denominator,
    'difference': str(r16-r1),
    'k1_source_sha256': hashlib.sha256((MAIN / 'rational_inward_candidate_results.json').read_bytes()).hexdigest(),
    'k16_source_sha256': hashlib.sha256((MAIN / 'rational_inward_k16_results.json').read_bytes()).hexdigest(),
    'self_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
assert result['k16_exact_rational_lower_strictly_exceeds_k1_exact_rational_lower']
assert result['k16_finite_L1_target_strictly_exceeds_k1_exact_rational_lower']
(HERE / 'compare_k16_to_k1_results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
