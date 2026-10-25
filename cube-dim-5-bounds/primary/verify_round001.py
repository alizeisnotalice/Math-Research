#!/usr/bin/env python3
"""Check integrated certificate identity and exact comparisons.
This does not replace the geometric reconstruction or analytic proofs.
"""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
r=Path(__file__).resolve().parents[1]
def read(p):return json.loads((r/p).read_text())
def sha(p):return hashlib.sha256((r/p).read_bytes()).hexdigest()
p=read('lower/rectangular_certificate.json')
a=read('audit/002_independent_lower_results.json')
l=read('results/bounds_ledger.json')
c=next(x for x in p['candidates'] if x['delta1']=='0' and x['delta2']=='1/64')
assert sha('lower/certify_rectangular.py')==p['source_sha256']
assert sha('audit/002_independent_lower_audit.py')==a['source_sha256']
assert c['coefficients']==a['ratio_coefficients']
assert sorted(c['retained_witness_counts'])==sorted(a['retained_witness_counts'])
assert c['union_engines_agree']
lo,hi=map(Q,c['ratio_interval']); alo,ahi=map(Q,a['ratio_interval'])
assert lo<hi and alo<ahi and lo<ahi and alo<hi
public=Q(read('literature/aldaz_a5_certificate.json')['interval'][1])
r4=Q(l['opening_thresholds']['lower']['exact_interval'][1])
assert min(lo,alo)>max(public,r4,Q(39863,20000))
assert l['lower']['exact_interval']==c['ratio_interval']
assert l['upper']['exact']=='32' and l['upper']['new_bound'] is False
print(json.dumps({'status':'PASS','scope':'identity, exact coefficients/witness equality, comparison only','retained_boxes':len(c['retained_witness_counts']),'strict_lower':'39863/20000','upper':'32','both_sides_breakthrough':False,'script_sha256':sha('primary/verify_round001.py')},indent=2))
