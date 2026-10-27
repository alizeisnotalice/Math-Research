#!/usr/bin/env python3
"""Portable exact rerun; requires only Python 3 standard library, writes nothing."""
from pathlib import Path
from fractions import Fraction
import hashlib,importlib.util,json
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
s=importlib.util.spec_from_file_location('cube8',H/'run.py');mod=importlib.util.module_from_spec(s);s.loader.exec_module(mod)
saved=json.loads((H/'results.json').read_text());reg=json.loads((H/'registration.json').read_text())
assert sha(H/'registration.json')==saved['registration_sha256']
assert sha(H/'run.py')==saved['code_sha256']
assert reg['ordered_m'][0]==64 and len(saved['trials'])==1
rr,tail,records,orbit=mod.build()
assert {str(k):str(v) for k,v in rr.items()}==saved['root_lower_numerators']
ref=Fraction(saved['comparison_upper']);inherited=Fraction(saved['inherited_finite_lower'])
actual=mod.calculate(64,rr,tail,records,orbit,ref,inherited); expected=saved['trials'][0]
for key in ['m','cell_upper','cell_prefilter_excluded','inner_states','outer_states','union_interval','states','positive_inward_raw_boxes','positive_outward_raw_boxes','containment_checks','minimum_width','finite_L1','winner']:
 assert actual[key]==expected[key],key
assert actual['winner'] and Fraction(actual['finite_L1']['margin'])>0
print('PASS: all roots, all labels/orbits, exact unions, containment, finite L1 and strict comparison')
