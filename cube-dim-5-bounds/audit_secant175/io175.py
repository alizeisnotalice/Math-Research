"""Own simple exact-file IO. No producer or mathematical evaluation imports."""
from pathlib import Path
import gzip,hashlib,json,os
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        while data:=f.read(1<<20):h.update(data)
    return h.hexdigest()
def read(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else Path(p).open()) as f:return json.load(f)
def lines(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else Path(p).open()) as f:
        for row in f:yield json.loads(row)
def put(p,obj):
    with Path(p).open('x') as f:json.dump(obj,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def lock(locks):
    for row in locks.values():assert sha(row['path'])==row['sha256'],row['path']
