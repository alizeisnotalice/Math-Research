import sys,json,multiprocessing as mp,time,hashlib
from pathlib import Path
sys.path.insert(0,'/Volumes/data/CodexProjects/T/Math-Research-strong-maximal/strong-maximal-endpoint/code')
from experiment92 import job
from verify85 import constants85
if __name__=='__main__':
 out=Path('/tmp/strong_fullcpu93');out.mkdir(exist_ok=True)
 mu=constants85()['mu_interval'];bs=('1/16','1/8','1/4','1/2','1','2','4','8','16','32','64','128')
 specs=[(k,b,mu) for k in range(64,4097,64) for b in bs];start=time.monotonic();pids=set();n=0
 with (out/'rows.jsonl').open('w') as f,mp.get_context('spawn').Pool(12) as pool:
  for r in pool.imap_unordered(job,specs,chunksize=1):
   f.write(json.dumps(r)+'\n');f.flush();n+=1;pids.add(r['pid'])
   if n%12==0:print(json.dumps({'completed':n,'total':len(specs),'workers':len(pids),'seconds':time.monotonic()-start}),flush=True)
 (out/'summary.json').write_text(json.dumps({'pass_batch':True,'cases':n,'actual_workers':len(pids),'pids':sorted(pids),'seconds':time.monotonic()-start,'scope':'Extended Gaussian near-one rank-window diagnostics only; no physical norm certificate.'},indent=2)+'\n')
