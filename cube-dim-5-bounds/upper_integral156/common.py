"""156 exact IO; no H/R/log/power evaluation."""
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

def dyadic(s):
    # Raw source fractions are already reduced; require a power-of-two denominator.
    parts=s.split('/');n=int(parts[0]);d=1 if len(parts)==1 else int(parts[1])
    assert len(parts)<=2 and d>0 and d&(d-1)==0
    return n,d.bit_length()-1

def endpoint_arrays(raw,times,point,key):
    assert raw['key']==key and raw['point']==point
    ends=raw['endpoints'];assert len(ends)==73 and [x['t'] for x in ends]==times
    pairs=[(dyadic(e['H']['lower']),dyadic(e['H']['upper'])) for e in ends]
    E=max(x[1] for pair in pairs for x in pair)
    lo=[a[0]<<(E-a[1]) for a,b in pairs];up=[b[0]<<(E-b[1]) for a,b in pairs]
    assert all(0<a<=b for a,b in zip(lo,up))
    assert all(lo[j]>(up[j+1] if j<72 else 0) for j in range(73))
    return lo,up,E

def lockcheck(locks):
    for v in locks.values():assert sha(v['path'])==v['sha256'],v['path']

def ceildiv(n,d):
    assert n>=0 and d>0
    return (n+d-1)//d
