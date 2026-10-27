#!/usr/bin/env python3
"""Check frozen saved results and exact comparisons; does not rerun union."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
lock=json.loads((H/'input_lock.json').read_text());reg=json.loads((H/'registration.json').read_text())
assert sha(H/'registration.json')==lock['registration_sha256']
assert sha(H/'geometry_r001.py')==lock['geometry_r001_sha256']
oldpath=H.parent/'round001/results.json';assert sha(oldpath)==lock['r001_results_sha256']
old=json.loads(oldpath.read_text());ref=F(old['trials'][0]['union_interval'][1]);inherit=F(old['trials'][0]['finite_L1']['lower']);rows=[]
for m in reg['candidates']:
 d=json.loads((H/f'candidate_{m}.json').read_text());c=d['candidate']
 assert d['m']==m and d['status']=='complete' and c['m']==m
 assert d['registration_sha256']==lock['registration_sha256'] and d['geometry_sha256']==lock['geometry_r001_sha256'] and d['r001_results_sha256']==lock['r001_results_sha256']
 assert F(d['comparison_upper'])==ref and F(d['inherited_lower'])==inherit
 assert c['states']<=reg['memo_limit_per_candidate'] and c['inner_states']+c['outer_states']==c['states']
 lo,hi=map(F,c['union_interval']);assert lo<=hi<=F(c['cell_upper']) and F(c['cell_upper'])>inherit
 winner=c.get('winner',False)
 if winner:
  finite=F(c['finite_L1']['lower']);assert lo>ref and finite>ref and F(c['finite_L1']['margin'])==finite-ref
 else: assert hi<ref
 rows.append({'m':m,'inner_vs_reference':(lo>ref)-(lo<ref),'outer_vs_reference':(hi>ref)-(hi<ref),'finite_winner':winner,'states':c['states']})
assert [r['m'] for r in rows if r['finite_winner']]==[96]
print(json.dumps({'status':'PASS saved arithmetic and hashes; no union rerun','results':rows},indent=2))
