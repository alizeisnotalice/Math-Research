#!/usr/bin/env python3
from pathlib import Path
import concurrent.futures,json,subprocess,sys,time,hashlib
H=Path(__file__).resolve().parent
def run(name):
 start=time.time();p=subprocess.run([sys.executable,'-B',str(H/'candidate.py'),name],capture_output=True,text=True);end=time.time()
 result={'name':name,'command':[sys.executable,'-B','candidate.py',name],'returncode':p.returncode,'started_unix':start,'ended_unix':end,'wall_seconds':end-start,'stdout':p.stdout,'stderr':p.stderr}
 with (H/f'receipt_{name}.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({'name':name,'returncode':p.returncode,'stdout':p.stdout.strip()}),flush=True)
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,['left','right']))
assert all(r['returncode']==0 for r in results)
