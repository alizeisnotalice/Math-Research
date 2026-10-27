import os,hashlib,json,gzip,time,sys
from pathlib import Path
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round129/recovered_continuation');D=Path(__file__).resolve().parent
start=time.time();count=total=empty=0
with gzip.open(D/'registered_inputs.jsonl.gz','xt',encoding='utf-8',compresslevel=6) as out:
 for root,dirs,files in os.walk(S,followlinks=False):
  dirs.sort();files.sort()
  for name in files:
   p=Path(root)/name
   if p.is_symlink() or not p.is_file():raise RuntimeError('nonregular '+str(p))
   st=p.stat(); count+=1;total+=st.st_size
   if count>200000 or total>2000000000:raise RuntimeError('cap exceeded')
   rel=p.relative_to(S).as_posix();h=hashlib.sha256()
   with p.open('rb') as f:
    while True:
     b=f.read(1024*1024)
     if not b:break
     h.update(b)
   now=p.stat()
   if now.st_size!=st.st_size or now.st_mtime_ns!=st.st_mtime_ns:raise RuntimeError('changed '+rel)
   e=rel.startswith('attempts_v129/') and rel.endswith('/stderr.txt')
   if e and st.st_size!=0:raise RuntimeError('nonempty stderr '+rel)
   empty+=e
   out.write(json.dumps({'path':rel,'bytes':st.st_size,'sha256':h.hexdigest(),'storage':'empty' if e else 'tar'},separators=(',',':'))+'\n')
(D/'registration_receipt.json').write_text(json.dumps({'status':'REGISTERED','files':count,'bytes':total,'empty_stderr':empty,'command':sys.argv,'exit_code':0,'elapsed_seconds':time.time()-start},indent=2)+'\n')
