"""No extraction or scientific recomputation; verify all byte-identical archive members."""
from pathlib import Path
import json,gzip,io,tarfile,hashlib
P=Path(__file__).resolve().parent
r=json.loads((P/'bundle_receipt.json').read_text())
for x in r['parts']:
 f=P/x['path'];assert f.stat().st_size==x['bytes']and hashlib.sha256(f.read_bytes()).hexdigest()==x['sha256']
index={x['archive_path']:x for x in map(json.loads,gzip.open(P/'inventory.jsonl.gz','rt'))if 'include_new_archive' in x['actions']}
class Reader(io.RawIOBase):
 def __init__(self):self.files=iter(r['parts']);self.f=None
 def readable(self):return True
 def read(self,n=-1):
  if n<0:n=1048576
  out=bytearray()
  while len(out)<n:
   if self.f is None:
    try:self.f=(P/next(self.files)['path']).open('rb')
    except StopIteration:break
   b=self.f.read(n-len(out));out.extend(b)
   if not b:self.f.close();self.f=None
  return bytes(out)
seen=set();total=0
with tarfile.open(fileobj=Reader(),mode='r|gz')as tf:
 for m in tf:
  assert m.isfile()and m.name in index and m.name not in seen
  b=tf.extractfile(m).read();x=index[m.name];assert len(b)==x['bytes']and hashlib.sha256(b).hexdigest()==x['sha256'];seen.add(m.name);total+=len(b)
assert len(seen)==r['files'] and total==r['source_bytes']
print('PASS',len(seen),'exact members; no scientific certification implied')
