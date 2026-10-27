from pathlib import Path
import gzip,json,hashlib,tarfile,io,os,time
p=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round127/sealed_bundle')
rows=[json.loads(x)for x in gzip.open(p/"inventory.jsonl.gz","rt")];selected=rows;index={x["archive_path"]:x for x in rows};assert len(index)==len(rows)
assert not list(p.glob("proof_chain*"))
start=time.monotonic()
class Parts:
 def __init__(self):self.f=None;self.n=0;self.size=0;self.records=[]
 def write(self,b):
  original=len(b)
  while b:
   if self.f is None:self.name=f'proof_chain.tar.gz.part{self.n:03d}';self.f=(p/self.name).open('xb');self.size=0;self.h=hashlib.sha256();self.n+=1
   part=b[:min(len(b),32*1024*1024-self.size)];self.f.write(part);self.h.update(part);self.size+=len(part);b=b[len(part):]
   if self.size==32*1024*1024:self.closepart()
  return original
 def flush(self):
  if self.f:self.f.flush()
 def closepart(self):
  if self.f:self.f.flush();os.fsync(self.f.fileno());self.f.close();self.records.append({'path':self.name,'bytes':self.size,'sha256':self.h.hexdigest()});self.f=None
sink=Parts()
with gzip.GzipFile(fileobj=sink,filename='',mode='wb',compresslevel=6,mtime=0)as gz:
 with tarfile.open(fileobj=gz,mode='w|',format=tarfile.PAX_FORMAT)as tf:
  for row in selected:
   data=Path(row['source_path']).read_bytes();assert len(data)==row['bytes']and hashlib.sha256(data).hexdigest()==row['sha256'],row['source_path']
   t=tarfile.TarInfo(row['archive_path']);t.size=len(data);t.mode=0o644;t.mtime=0;tf.addfile(t,io.BytesIO(data))
sink.closepart()
class ReadParts(io.RawIOBase):
 def __init__(self):self.files=iter(sink.records);self.f=None
 def readable(self):return True
 def read(self,n=-1):
  if n<0:n=1024*1024
  out=bytearray()
  while len(out)<n:
   if self.f is None:
    try:self.f=(p/next(self.files)['path']).open('rb')
    except StopIteration:break
   b=self.f.read(n-len(out));out.extend(b)
   if not b:self.f.close();self.f=None
  return bytes(out)
verified=set();total=0
with tarfile.open(fileobj=ReadParts(),mode='r|gz')as tf:
 for item in tf:
  assert item.isfile()and item.name in index and item.name not in verified;data=tf.extractfile(item).read();row=index[item.name];assert len(data)==row['bytes']and hashlib.sha256(data).hexdigest()==row['sha256'];verified.add(item.name);total+=len(data)
assert len(verified)==len(index) and total==sum(x['bytes']for x in selected)
(p/'bundle_receipt.json').write_text(json.dumps({'status':'BUNDLE_AND_FULL_STREAM_HASH_CHECK_PASS','files':len(verified),'source_bytes':total,'parts':sink.records,'wall_seconds':time.monotonic()-start,'scientific_calls':0,'source_mutations':0,'excludes':'all selected C127 crash/recovery replay records preserved; prior12081/runtime references not duplicated; live129 excluded'},indent=2)+'\n')
print(json.dumps({'status':'PASS','files':len(verified),'source_bytes':total,'compressed_bytes':sum(x['bytes']for x in sink.records),'parts':len(sink.records),'wall_seconds':time.monotonic()-start}))
