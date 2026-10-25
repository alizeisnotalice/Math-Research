from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
r=Path(__file__).resolve().parents[2];d=r/'lower/round003'
paths=[d/'second_spacing_certificate_results.json',d/'second_spacing_k16_certificate_results.json',r/'audit/lower003/rational_inward_k16_results.json']
a,b,c=[json.loads(p.read_text()) for p in paths]
assert Q(c['exact_rational_lower'])>Q(a['exact_interval'][1])
assert Q(b['exact_interval'][0])>Q(a['exact_interval'][1])
out={'status':'PASS_exact_current_k1_comparison','k16_inward_lower':c['exact_rational_lower'],'current_k1_upper':a['exact_interval'][1],'k16_full_lower':b['exact_interval'][0],'inputs':{str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(out['status'])
