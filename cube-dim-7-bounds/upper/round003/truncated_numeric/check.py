#!/usr/bin/env python3
"""Portable integrity checker; stdlib only, no scientific kernel recomputation."""
import json, hashlib, sys
from pathlib import Path
p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
m=json.loads((p/'manifest.json').read_text())
for name,digest in m['files'].items():
 if sha(p/name)!=digest:raise SystemExit('HASH MISMATCH: '+name)
r=json.loads((p/'results.json').read_text()); tasks=[json.loads(x) for x in (p/'task-ledger.jsonl').read_text().splitlines()]
assert r['status']=='completed' and len(tasks)==257 and sorted(x['id'] for x in tasks)==list(range(257))
assert all(x['status']=='completed' for x in tasks) and sum(x['points'] for x in tasks)==4104
assert sum(x['point_epsilon_shape_terms'] for x in tasks)==4104*4*8192==r['task_point_epsilon_shape_terms']
assert all(x['min_q']>0 for x in tasks) and tasks[0]['warm_cholesky_max_relative_error']<1e-8
assert len(r['minima'])==4 and len(r['old_S3_points'])==8 and sha(p/'points.npy')==r['point_file_sha256']
assert all(len(x['block_means'])==8 and len(x['point_floathex'])==7 for x in r['minima'])
if len(sys.argv)>1:
 assert sha(Path(sys.argv[1]))==r['cache_sha256']; print('Optional cache SHA256 verified.')
print('PASS: manifest, unique terminal tasks, registered counts, Q diagnostics, and retained results. This is not a proof or a recomputation of H.')
