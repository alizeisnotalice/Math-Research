#!/usr/bin/env python3
"""Check portable artifacts and conservation identities, not sampling accuracy."""
import json,hashlib,math
from pathlib import Path
from fractions import Fraction
p=Path(__file__).resolve().parent;sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();m=json.loads((p/'manifest.json').read_text())
for f,h in m['files'].items():assert sha(p/f)==h, f
pts=json.loads((p/'points.json').read_text())['points'];x=[Fraction(z) for z in pts[0]['exact_rational']]; near=[Fraction(z) for z in pts[1]['exact_rational']];far=[Fraction.from_float(float.fromhex(z)) for z in pts[2]['exact_floathex']]
assert x==[1,1,1,1,0,1,1] and near==[Fraction(1)]+[(1-Fraction(7,1000003))*x[j]+Fraction(j,1000003) for j in range(1,7)]
assert len({abs(z) for z in near})==7 and len({abs(z) for z in far})==7
r=json.loads((p/'results.json').read_text());ts=[json.loads(z) for z in (p/'terminal.jsonl').read_text().splitlines()];assert r['status']=='completed' and len(ts)==15 and {x['handle'] for x in ts}=={f'p{i}s{j}' for i in range(3) for j in range(5)}
for a in ts:
 i,j=int(a['handle'][1]),int(a['handle'][3]);assert a['status']=='completed' and a['seed']==[2026102704,i,j] and a['n']==2**20 and a['chunks']==1024 and len(a['block_means'])==16
 assert all(sum(row)==a['n'] for row in a['kernel_bin_counts']) and sum(a['u_bin_counts'])==a['n']
 for k in range(7):assert math.isclose(sum(v[k] for v in a['u_bin_sums']),a['sum'][k],rel_tol=3e-12,abs_tol=1e-18)
 for k in range(4):assert math.isclose(sum(a['kernel_bin_sums'][k]),a['sum'][k],rel_tol=3e-12,abs_tol=1e-18)
 for k in range(3):assert math.isclose(a['sum'][k+1]-a['sum'][k],a['sum'][k+4],rel_tol=1e-10,abs_tol=1e-10)
print('PASS: hashes, exact point membership tests, 15 terminal seeds and counts, histogram and time-band conservation. No expectation error is certified.')
