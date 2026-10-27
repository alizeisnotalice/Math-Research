import os, sys, json, hashlib, tarfile, gzip, time, subprocess
from pathlib import Path
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round129/recovered_continuation')
D=Path(__file__).resolve().parent
LIMIT=32*1024*1024
class Split:
 def __init__(self): self.i=0;self.n=0;self.f=None;self.parts=[]
 def write(self,b):
  total=len(b)
  while b:
   if self.f is None:
    p=D/f'proof_chain.tar.gz.part{self.i:03}';self.f=p.open('xb');self.h=hashlib.sha256();self.n=0
   chunk=b[:LIMIT-self.n];self.f.write(chunk);self.h.update(chunk);self.n+=len(chunk);b=b[len(chunk):]
   if self.n==LIMIT:self.finish()
  return total
 def finish(self):
  if self.f:self.f.close();self.parts.append({'path':f'proof_chain.tar.gz.part{self.i:03}','bytes':self.n,'sha256':self.h.hexdigest()});self.f=None;self.i+=1
 def flush(self):
  if self.f:self.f.flush()
class Reader:
 def __init__(self,f):self.f=f;self.h=hashlib.sha256();self.n=0
 def read(self,n=-1):
  b=self.f.read(n);self.h.update(b);self.n+=len(b);return b
start=time.time(); paths=[];size=0
with gzip.open(D/'registered_inputs.jsonl.gz','rt',encoding='utf8') as f: registered={r['path']:r for r in map(json.loads,f)}
for root,dirs,files in os.walk(S,followlinks=False):
 dirs.sort();files.sort()
 for name in files:
  p=Path(root)/name
  if p.is_symlink() or not p.is_file():raise RuntimeError('nonregular '+str(p))
  st=p.stat();paths.append((p,st.st_size,st.st_mtime_ns));size+=st.st_size
  if len(paths)>200000 or size>2000000000:raise RuntimeError('input cap exceeded')
split=Split()
with gzip.open(D/'inventory.jsonl.gz','xt',encoding='utf8',compresslevel=6) as inv:
 with gzip.GzipFile(fileobj=split,mode='wb',compresslevel=6,mtime=0) as gz:
  with tarfile.open(fileobj=gz,mode='w|',format=tarfile.PAX_FORMAT) as tar:
   for p,n,mt in paths:
    rel=p.relative_to(S).as_posix()
    if rel not in registered: raise RuntimeError('unregistered source '+rel)
    expected=registered.pop(rel)
    if expected['bytes'] != n: raise RuntimeError('registered size mismatch '+rel)
    if rel.startswith('attempts_v129/') and rel.endswith('/stderr.txt'):
     if n != 0 or expected['sha256'] != hashlib.sha256(b'').hexdigest() or expected['storage'] != 'empty': raise RuntimeError('stderr registration mismatch '+rel)
     with p.open('rb') as f:
      if f.read(1): raise RuntimeError('stderr changed '+rel)
     st=p.stat()
     if st.st_size != n or st.st_mtime_ns != mt: raise RuntimeError('source changed '+rel)
     inv.write(json.dumps({'path':rel,'bytes':0,'sha256':hashlib.sha256(b'').hexdigest(),'storage':'empty'},separators=(',',':'))+'\n')
     continue
    info=tarfile.TarInfo(rel);info.size=n;info.mode=0o644;info.mtime=0
    with p.open('rb') as f:
     reader=Reader(f);tar.addfile(info,reader)
    st=p.stat()
    if reader.n!=n or st.st_size!=n or st.st_mtime_ns!=mt:raise RuntimeError('source changed '+rel)
    if expected['sha256'] != reader.h.hexdigest() or expected['storage'] != 'tar': raise RuntimeError('registered hash mismatch '+rel)
    inv.write(json.dumps({'path':rel,'bytes':n,'sha256':reader.h.hexdigest()},separators=(',',':'))+'\n')
split.finish()
if registered: raise RuntimeError('registered source missing')
(D/'build_receipt.json').write_text(json.dumps({'status':'BUILT_NOT_YET_VERIFIED','files':len(paths),'source_bytes':size,'parts':split.parts,'elapsed_seconds':time.time()-start,'command':sys.argv,'exit_code':0},indent=2)+'\n')
