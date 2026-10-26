#!/usr/bin/env python3
"""Explicitly rebuild one080 cell in a NEW directory, reusing old066 endpoint."""
from pathlib import Path
import argparse,gzip,json,importlib.util,hashlib,time
H=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--rank',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert 0<=a.rank<52360;assert not a.output.exists(),'new output directory required; never overwrite evidence';a.output.mkdir(parents=True)
spec=importlib.util.spec_from_file_location('engine080',H/'engine_v2.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.S=a.output.resolve();reg=json.loads((H/'registration.json').read_text());policy=json.loads((H/'execution_policy_v2.json').read_text());fp=policy['old_fingerprint'] if a.rank<policy['transition_rank'] else hashlib.sha256((m.sha(H/'engine_v2.py')+m.sha(H/'registration.json')+m.sha(H/'execution_policy_v2.json')).encode()).hexdigest();table=json.loads(gzip.decompress((H/'geometry017.json.gz').read_bytes()));root=json.loads(gzip.decompress((H/'roots.json.gz').read_bytes()))[a.rank];m.init_worker(table,reg['time_endpoints'],fp);result=m.worker(root);print(json.dumps(result,indent=2))
if result['status']=='complete' and (H/'initial_cells.json.gz').exists():
 compact=json.loads(gzip.decompress((H/'initial_cells.json.gz').read_bytes()));expected=dict(zip(compact['columns'],compact['rows'][a.rank]))['mathematical_sha256'];ev=json.loads(gzip.decompress(m.evidencepath(a.rank).read_bytes()));ev.pop('runtime');actual=hashlib.sha256(json.dumps(ev,sort_keys=True,separators=(',',':')).encode()).hexdigest();assert actual==expected;print('deterministic mathematical payload SHA matches',actual)
