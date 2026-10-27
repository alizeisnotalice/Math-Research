"""Portable175 byte/JSON IO; no mathematical consumer imported."""
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

def lines(p):
    with gzip.open(p,'rt') as f:
        for line in f:yield json.loads(line)

def put(p,v):
    with Path(p).open('x') as f:json.dump(v,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())

def ref(p):
    p=Path(p);return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}

def checked(r):
    p=Path(r['path']);assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'];return p

def gzwrite(p,raw):
    with Path(p).open('xb') as f:
        with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0,compresslevel=6) as g:g.write(raw)

def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':')).encode('utf-8')
