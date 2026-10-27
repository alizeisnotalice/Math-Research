#!/usr/bin/env python3
"""Exact-dyadic 5040-chunk coverage/reuse and outward sum checker; no phase re-evaluation."""
from pathlib import Path
from fractions import Fraction
import json,hashlib
p=Path(__file__).resolve().parent;sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();m=json.loads((p/'manifest.json').read_text())
for n,h in m['files'].items():assert sha(p/n)==h,n
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,n=line.split('  ',1);q=Path(n);assert sha(q if q.is_absolute() else p/q)==h,n
oldmanifest=json.loads((p/'r011-manifest.json').read_text());mp={'kernel.py':'run.py','r011-results.json':'results.json','cells.jsonl':'cells.jsonl','r009-rays.json':'r009-rays.json','r009-inputs.json':'r009-inputs.json','r010-manifest.json':'r010-manifest.json','testcases.json':'testcases.json'}
for a,b in mp.items():assert sha(p/a)==oldmanifest['files'][b]
pilot=json.loads((p/'r011-results.json').read_text());assert pilot['status']=='completed' and pilot['generic_planning_width_pass'] and not pilot['promotions'];old={x['chunk']:x for x in pilot['generic_chunks']};assert set(old)==set(range(4))
R=json.loads((p/'result.json').read_text());receipt=json.loads((p/'receipt.json').read_text());T=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];assert len(T)==5040 and {x['chunk'] for x in T}==set(range(5040)) and all(x['status']=='completed' for x in T)
assert [x['chunk'] for x in T[:4]]==[0,1,2,3] and receipt['status']==R['status']=='completed'
def val(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
low=Fraction();up=Fraction();sums={k:0 for k in R['counts']}
for x in T:
 j=x['chunk'];assert x['precision']==256 and x['start_phase']==64*j and x['end_phase']==64*(j+1)
 assert x['source']==('reused_r011' if j<4 else 'new')
 if j<4:assert x['interval']==old[j]['interval']
 lo=val(x['interval']['lower']);hi=val(x['interval']['upper']);assert lo<=hi;low+=lo;up+=hi
 for k in sums:sums[k]+=x['counts'][k]
 assert x['counts']['top_requests']==64*212
assert sums==R['counts'] and sums['top_requests']==68382720 and R['halfphases']==322560 and R['new_chunks']==5036 and R['reused_chunks']==4
L=val(R['full_H']['lower']);U=val(R['full_H']['upper']);assert 0<L<=low<=up<=U and U-L<=Fraction(1,10**12) and val(R['width']['upper'])>=U-L
assert R['point_floathex']==json.loads((p/'testcases.json').read_text())['generic_case']['x_floathex'] and R['epsilon']=='1/1024' and R['positive'] and R['width_pass'] and R['not_a_global_g_certificate']
print('PASS: immutable pilot source and 4 unchanged reused chunks, 5036 new chunk terminals, exact 322560 phase coverage/counts and outward exact-dyadic full sum. No additional phases evaluated; single point does not give g.')
