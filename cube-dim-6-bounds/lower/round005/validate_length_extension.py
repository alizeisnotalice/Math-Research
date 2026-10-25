"""Bind exact span, corrected direct-containment and endpoint-grid certificates."""
from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parents[2]
span=r/'audit/lower005/span_prune_results.json';direct=r/'audit/lower005/containment_check_corrected_results.json';grid=r/'lower/round005/endpoint_grid_cover_fixed_results.json'
s,d,g=[json.loads(p.read_text()) for p in [span,direct,grid]]
assert s['source_sha256']==hashlib.sha256((span.parent/'span_prune.py').read_bytes()).hexdigest()
assert d['source_sha256']==hashlib.sha256((direct.parent/'containment_check.py').read_bytes()).hexdigest()
assert g['code_sha256']==hashlib.sha256((grid.parent/'endpoint_grid_cover_fixed.py').read_bytes()).hexdigest()
assert d['status']=='strict_interval_containment_audit_corrected_radius' and g['status']=='pass'
assert s['m5']['new_shell_total']==5**6-4**6==11529
assert s['m6']['new_shell_total']==6**6-5**6==31031
for key in ['m5','m6']:
 assert {tuple(e) for e in s[key]['new_shell_survivors']}=={tuple(x['ell']) for x in d['shells'][key]}
 assert all(x['new_canonical_half_box_certified_nonempty'] for x in d['shells'][key])
missing={tuple(x['ell']) for x in d['shells']['m5'] if not x['contained_in_old_lte4']}
assert len(missing)==5 and missing=={tuple(x['ell']) for x in g['targets']}
assert all(x['status']=='covered' for x in g['targets'])
assert len(d['shells']['m5'])==17 and len(d['shells']['m6'])==7
assert all(x['contained_in_old_lte4'] for x in d['shells']['m6'])
out={'status':'PASS_strict_finite_canonical_length_extension_no_gain','cap5_new_vectors':11529,'cap5_positive_survivors':17,'cap5_directly_contained':12,'cap5_union_contained':5,'cap6_new_vectors':31031,'cap6_positive_survivors_all_directly_contained':7,'conclusion':'For the fixed R6-003 spacing and canonical half-cell witness family, extending ell cap4 to5 or6 leaves the positive-volume union unchanged. No statement for ell>6, shifted witnesses, other lattices, weights, or C6 optimality.','inputs':{str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [span,direct,grid]},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
