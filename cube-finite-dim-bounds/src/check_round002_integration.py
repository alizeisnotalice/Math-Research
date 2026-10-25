"""Check packaging identities and strict comparisons; does not rerun search or union geometry."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,ast
p=Path(__file__).resolve().parents[1];root=p.parent
base=json.loads((p/'baseline/original_certificates.json').read_text())
n3=json.loads((p/'lower/n3/certificate_result.json').read_text())
assert F(n3['ratio_interval'][1])<F(base[0]['lower_bound_interval'][0])
assert n3['b']==2 and n3['critical_inner_radius_count']==2
assert hashlib.sha256((p/'lower/n3/certificate_input.json').read_bytes()).hexdigest()==n3['input_sha256']
n4=json.loads((p/'lower/n4/scan_certificate.json').read_text())
for row in n4['rows']+n4['enriched_rows']:
 assert F(row['outer_ratio_upper'])<F(base[1]['lower_bound_interval'][0])
checks=0
m=json.loads((p/'lower/n3/source_manifest.json').read_text())
for r in m['entries']:
 if r['path'].startswith('cube-finite-dim-bounds/lower/n3/'):
  assert hashlib.sha256((root/r['path']).read_bytes()).hexdigest()==r['sha256'];checks+=1
m=json.loads((p/'lower/n4/source_manifest.json').read_text())
for name,r in m['files'].items():
 assert hashlib.sha256((p/'lower/n4'/name).read_bytes()).hexdigest()==r['sha256'];checks+=1
for d in ['src','lower/n3','lower/n4','upper','audit']:
 for f in (p/d).glob('*.py'):ast.parse(f.read_text())
print(json.dumps({'status':'PASS','own_file_hashes_checked':checks,'n3_exact_comparisons':1,'n4_exact_comparisons':len(n4['rows'])+len(n4['enriched_rows']),'scope':'packaging, syntax, stored exact comparisons; geometry audited separately'},indent=2))
