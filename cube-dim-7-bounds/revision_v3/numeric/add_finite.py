import sys
if hasattr(sys,"set_int_max_str_digits"): sys.set_int_max_str_digits(0)
"""Add one finite-a certificate to each previously selected rational t."""
import csv,json
from pathlib import Path
from fractions import Fraction as Q
from verify import finite_certify
p=Path(__file__).parent
rows=json.loads((p/'table.json').read_text())
for r in rows:
    vals=finite_certify(r['n'],Q(r['t']),Q(r['B_upper']))
    assert vals[-1]<Q(r['U'])
    r.update(zip(('finite_A','finite_delta_upper','finite_E_upper','finite_tau','finite_bound'),map(str,vals)))
    r['finite_a']='2^256'
(p/'table.json').write_text(json.dumps(rows,indent=2)+'\n')
with (p/'table.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print('PASS: finite-a rational bounds below U for',len(rows),'dimensions')
