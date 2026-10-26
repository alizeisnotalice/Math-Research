from pathlib import Path
import hashlib,json,os
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def atomic(p,x,fresh=False):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp')
 if fresh:assert not p.exists()
 with open(tmp,'w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
