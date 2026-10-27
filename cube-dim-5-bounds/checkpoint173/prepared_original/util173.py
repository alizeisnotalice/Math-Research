from pathlib import Path
import hashlib,json,os,gzip
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read(p):
 with gzip.open(p,'rt') if str(p).endswith('.gz') else open(p) as f:return json.load(f)
def atomic(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
 with open(tmp,'w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def create(p,x):assert not Path(p).exists();atomic(p,x)
def append(f,row):
 row={**row,'hash':digest(row)};f.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode());f.flush();os.fsync(f.fileno());return row['hash']
