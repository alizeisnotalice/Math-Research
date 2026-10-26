#!/usr/bin/env python3
"""Explicit rebuild of one registered082 child, in a fresh output directory."""
from pathlib import Path
import argparse,gzip,json,hashlib
import kernel as K
H=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--key',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists();a.output.mkdir(parents=True);K.S=a.output.resolve();reg=json.loads((H/'registration.json').read_text());fp=hashlib.sha256((K.sha(H/'manager.py')+K.sha(H/'kernel.py')+K.sha(H/'registration.json')).encode()).hexdigest();jobs=[]
for path in sorted(H.glob('batch*_jobs.json.gz')):jobs+=json.loads(gzip.decompress(path.read_bytes()))
row=next(j for j in jobs if j['key']==a.key);table=json.loads(gzip.decompress((H.parent/'upper_global080/geometry017.json.gz').read_bytes()));K.init_worker(table,reg['time_endpoints'],fp);result=K.worker(row);print(json.dumps(result,indent=2));assert result['status']=='complete';ev=json.loads(gzip.decompress(K.evidencepath(a.key).read_bytes()));ev.pop('runtime');expected=json.loads(gzip.decompress((H/'mathematical_hashes.json.gz').read_bytes()))[a.key];assert hashlib.sha256(K.serial(ev).encode()).hexdigest()==expected;print('deterministic mathematical SHA matches')
