#!/usr/bin/env python3
from pathlib import Path
import concurrent.futures,json,subprocess,sys,time
H=Path(__file__).resolve().parent
reg=json.loads((H/'registration.json').read_text())
def run(event):
 rank=event['rank'];start=time.time();p=subprocess.run([sys.executable,'-B',str(H/'candidate.py'),str(rank)],capture_output=True,text=True);end=time.time()
 receipt={'rank':rank,'command':[sys.executable,'-B','candidate.py',str(rank)],'returncode':p.returncode,'started_unix':start,'ended_unix':end,'wall_seconds':end-start,'stdout':p.stdout,'stderr':p.stderr,'peak_rss_bytes':None}
 try:receipt['peak_rss_bytes']=json.loads((H/f'candidate_{rank:02}.json').read_text())['peak_rss_bytes']
 except Exception as exc:receipt['rss_unavailable_reason']=repr(exc)
 with (H/f'receipt_{rank:02}.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
 print(json.dumps({'rank':rank,'returncode':p.returncode,'stdout':p.stdout.strip(),'peak_rss_bytes':receipt['peak_rss_bytes']}),flush=True)
 return receipt
with concurrent.futures.ThreadPoolExecutor(max_workers=reg['parallel_limit']) as pool:results=list(pool.map(run,reg['selected']))
with (H/'runner_status.json').open('x') as f:json.dump({'all_children_returned_zero':all(r['returncode']==0 for r in results),'children':len(results)},f,indent=2);f.write('\n')
