#!/usr/bin/env python3
"""Integrity and bookkeeping only; no Monte Carlo replay or expectation certificate."""
from pathlib import Path
import json,hashlib,struct,math,gzip
p=Path(__file__).resolve().parent;sha=lambda n:hashlib.sha256((p/n).read_bytes()).hexdigest();manifest=json.loads((p/'manifest.json').read_text())
for n,h in manifest['files'].items():assert sha(n)==h,n
assert hashlib.sha256(gzip.decompress((p/'terminal.jsonl.gz').read_bytes())).hexdigest()==json.loads((p/'source-manifest.json').read_text())['files']['terminal.jsonl']
r=json.loads((p/'result_metadata.json').read_text());T=[json.loads(x) for x in gzip.decompress((p/'terminal.jsonl.gz').read_bytes()).decode().splitlines()];S=json.loads((p/'summary.json').read_text());assert r['status']=='completed' and len(T)==15 and {x['handle'] for x in T}=={f's{i:02}' for i in range(15)}
assert all(x['status']=='completed' and x['completed_chunks']==1024 and x['completed_shapes']==2**20 and len(x['subblock_means'])==8 for x in T)
assert all(x['seed']==[2026102707,int(x['handle'][1:])] for x in T) and sum(x['completed_shapes'] for x in T)==15728640
fallbacks=[f for x in T for f in x['fallbacks']];assert len(fallbacks)==r['fallback_count']==S['fallback_count']<=4096
for f in fallbacks:
 assert f['status']=='completed' and len(f['gamma_floathex'])==64 and 0<=f['chunk']<1024 and 0<=f['offset']<1024
 values=[float.fromhex(v) for v in f['gamma_floathex']];assert all(math.isfinite(v) and v>0 for v in values)
 assert hashlib.sha256(struct.pack('<64d',*values)).hexdigest()==f['gamma_sha256']
assert [sum(t['raw_violations'][i] for t in T) for i in range(8)]==S['raw_violations']
assert [sum(t['tolerance_violations'][i] for t in T) for i in range(8)]==S['tolerance_violations']
for j in range(21):
 total=sum(x['sum'][j] for x in T);mean=total/15728640;assert math.isclose(mean,S['stats'][j]['mean'],rel_tol=2e-15,abs_tol=1e-18)
 for t in T:assert math.isclose(sum(b[j] for b in t['subblock_means'])/8,t['sum'][j]/t['completed_shapes'],rel_tol=2e-13,abs_tol=1e-17)
for i in range(6):assert 0<=S['stats'][i]['minimum']<=S['stats'][i]['maximum']<=r['upper_envelope']*(1+1e-14)
assert S['stats'][13]['minimum']<0<S['stats'][13]['maximum']
print('PASS: hashes, 15 terminal streams/counts, 1997 preserved same-input fallback hashes, paired summaries and block arithmetic; empirical diagnostics only, no expectation or I certificate.')
