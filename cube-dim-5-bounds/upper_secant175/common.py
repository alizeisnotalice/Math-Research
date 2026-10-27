"""175 exact metadata IO only; never imports an older mathematical consumer."""
from pathlib import Path
import gzip,hashlib,json,os

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)

def put(p,x,atomic=False):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp') if atomic else p
    with tmp.open('w' if atomic else 'x') as f:
        json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    if atomic:os.replace(tmp,p)

def lines(p):
    with gzip.open(p,'rt') as f:
        for line in f:yield json.loads(line)

def lockcheck(locks):
    for x in locks.values():
        p=Path(x['path']);assert p.stat().st_size==x['bytes'] and sha(p)==x['sha256'],str(p)

def inputcheck(src):
    lockcheck(src['source_locks'])
    pre=read(src['source_locks']['metadata166']['path'])
    for x in pre['shards']+[pre['topology']]:assert sha(x['path'])==x['sha256']
    root=Path(src['source_locks']['completion166']['path']).parent
    for name,x in read(src['source_locks']['outputs166']['path'])['files'].items():
        p=root/name;assert p.stat().st_size==x['bytes'] and sha(p)==x['sha256']
