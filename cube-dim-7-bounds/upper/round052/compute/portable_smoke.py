#!/usr/bin/env python3
"""Read-only copy of the whitelist to a clean temporary tree and run portable check."""
import os,sys,tempfile,subprocess,shutil,json,time,hashlib
from pathlib import Path
D=Path(__file__).resolve().parent;start=time.monotonic()
names=['R052_REPORT.md','R052_REPORT_producer_original.md','ACTIVATION.md','ARCHIVE.md','registration.json','results.json','r052_integer.py','equivalence.json','verify_top64.py','verify_top64.json','monitor.jsonl','portable_check.py','portable_smoke.py','portable_smoke_result.json','archive_checksums.sha256']
with tempfile.TemporaryDirectory(prefix='r052-portable-') as t:
 dest=Path(t)/'cube-dim-7-bounds/upper/round052/compute';dest.mkdir(parents=True)
 for name in names:shutil.copy2(D/name,dest/name)
 for line in (dest/'archive_checksums.sha256').read_text().splitlines():
  digest,name=line.split(None,1);assert hashlib.sha256((dest/name.strip()).read_bytes()).hexdigest()==digest
 p=subprocess.run([sys.executable,str(dest/'portable_check.py')],capture_output=True,text=True,check=True)
 result=json.loads(p.stdout);assert result['status']=='PASS'
 print(json.dumps({'status':'PASS','portable_checker':result,'copied_files':len(names),'wall_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'temporary read-only source copy; no grid enumeration'},indent=2))
