import sys,json,hashlib,tarfile,gzip,time
from pathlib import Path
D=Path(__file__).resolve().parent
start=time.time();r=json.loads((D/'build_receipt.json').read_text())
class Parts:
 def __init__(self,parts):self.parts=iter(parts);self.f=None;self.h=None;self.n=0;self.part=None
 def read(self,n=-1):
  out=[]; remain=n
  while remain!=0:
   if self.f is None:
    try:self.part=next(self.parts)
    except StopIteration:break
    self.f=(D/self.part['path']).open('rb');self.h=hashlib.sha256();self.n=0
   b=self.f.read(remain)
   if b:self.h.update(b);self.n+=len(b);out.append(b);remain-=len(b)
   else:
    self.f.close();self.f=None
    if self.n!=self.part['bytes'] or self.h.hexdigest()!=self.part['sha256']:raise RuntimeError('part mismatch')
  return b''.join(out)
p=Parts(r['parts']); count=0;size=0
with gzip.open(D/'inventory.jsonl.gz','rt',encoding='utf8') as inv:
 with gzip.GzipFile(fileobj=p,mode='rb') as gz:
  with tarfile.open(fileobj=gz,mode='r|') as tar:
   for member in tar:
    row=json.loads(next(inv))
    while row.get('storage') == 'empty':
     if row['bytes'] != 0 or row['sha256'] != hashlib.sha256(b'').hexdigest() or not row['path'].startswith('attempts_v129/') or not row['path'].endswith('/stderr.txt'): raise RuntimeError('invalid empty entry')
     count+=1
     row=json.loads(next(inv))
    h=hashlib.sha256();f=tar.extractfile(member)
    while True:
     b=f.read(1024*1024)
     if not b:break
     h.update(b)
    if member.name!=row['path'] or member.size!=row['bytes'] or h.hexdigest()!=row['sha256']:raise RuntimeError('member mismatch '+member.name)
    count+=1;size+=member.size
   for line in inv:
    row=json.loads(line)
    if row.get('storage') != 'empty' or row['bytes'] != 0 or row['sha256'] != hashlib.sha256(b'').hexdigest() or not row['path'].startswith('attempts_v129/') or not row['path'].endswith('/stderr.txt'): raise RuntimeError('unmatched inventory')
    count+=1
  while gz.read(1024*1024):pass
while p.read(1024*1024):pass
if count!=r['files'] or size!=r['source_bytes']:raise RuntimeError('count mismatch')
(D/'verify_receipt.json').write_text(json.dumps({'status':'FULL_STREAM_SHA256_PASS','files':count,'source_bytes':size,'elapsed_seconds':time.time()-start,'command':sys.argv,'exit_code':0},indent=2)+'\n')
