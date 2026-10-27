#!/usr/bin/env python3
"""Read-only integrity and saved-result consistency, not a scientific rerun."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,sys
D=Path(__file__).resolve().parent
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 m=json.loads((D/'manifest.json').read_text())
 assert all(h(D/n)==v for n,v in m['files'].items())
 lock=json.loads((D/'sourcehash.json').read_text())
 assert h(D/'run.py')==lock['run_sha256']
 for root,files in [('theory_dir','theory'),('r010_dir','r010')]:
  for n,v in lock[files].items():assert h(Path(lock[root])/n)==v
 r=json.loads((D/'result.json').read_text());assert r['status']=='PASS' and r['orbit_id']=='psi7-004'
 assert r['det_L']==-11008 and r['incident_ordered_simplices']==36 and r['stabilizer_size']==8
 assert len(r['gates'])==3
 for g in r['gates']:
  assert Fraction(g['template_value'])==Fraction(g['star_value'])>0
  assert Fraction(g['x_dot_v'])==1 and Fraction(g['min_star_denominator'])>0
  assert g['star_terms']==288 and g['template_paired_terms']==10 and g['all_denominators_positive'] and g['exact_equal']
 print(json.dumps({'status':'PASS','scope':'hash and saved data consistency only; no scientific recomputation'}))
if __name__=='__main__':
 try:main()
 except Exception as e:print(json.dumps({'status':'FAIL','error':repr(e)}));sys.exit(1)
