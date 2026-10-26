"""Protocol helpers; verified import of unchanged131 exact corner core."""
from pathlib import Path
import hashlib,importlib.util,json,os
CORE_PATH=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_crossanchor131/consume.py')
CORE_SHA='071380d9de8a8376f9e92c8b6f2ba40c8fe3c9d5ed73d96b00f78fad294f6eb8'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def write(p,x):
 p=Path(p);t=p.with_suffix(p.suffix+'.tmp')
 with open(t,'w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def load_core():
 assert sha(CORE_PATH)==CORE_SHA
 spec=importlib.util.spec_from_file_location('verified131_exact_core',CORE_PATH);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 assert sha(CORE_PATH)==CORE_SHA;return m.single_source_certificate
