"""Single owned scientific child; hard180 seconds includes kill/reap."""
import json,subprocess,time,os,sys
from pathlib import Path
H=Path(__file__).resolve().parent;target=H/'process_receipt.json'
if target.exists():raise RuntimeError('Single attempt only')
start=time.monotonic()
p=subprocess.Popen([sys.executable,'-B',str(H/'lower_window.py')],cwd=H,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
r={'guard_pid':os.getpid(),'owned_pid':p.pid,'status':'RUNNING','wall_limit_including_cleanup':180,'automatic_retries':0}
target.write_text(json.dumps(r,indent=2)+'\n')
try:out,_=p.communicate(timeout=178)
except subprocess.TimeoutExpired:
    p.kill();out,_=p.communicate(timeout=2);r['timeout']=True
r.update(status='EXITED_AND_REAPED',exit_code=p.returncode,total_wall_seconds=time.monotonic()-start)
target.write_text(json.dumps(r,indent=2)+'\n')
print(out,end='');print(json.dumps(r));raise SystemExit(p.returncode)
