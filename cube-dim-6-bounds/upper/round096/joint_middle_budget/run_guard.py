"""One owned child,180 seconds including forced termination and reaping."""
import json,subprocess,time,sys,os
from pathlib import Path
H=Path(__file__).resolve().parent;R=H/'process_receipt.json'
if R.exists():raise RuntimeError('One invocation only')
start=time.monotonic();p=subprocess.Popen([sys.executable,'-B',str(H/'evaluate.py')],cwd=H,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
r={'status':'RUNNING','owned_pid':p.pid,'guard_pid':os.getpid(),'total_wall_limit_seconds':180,'automatic_retries':0}
R.write_text(json.dumps(r,indent=2)+'\n')
try:out,_=p.communicate(timeout=178)
except subprocess.TimeoutExpired:
    p.kill();out,_=p.communicate(timeout=2);r['timeout']=True
r.update(status='EXITED_AND_REAPED',exit_code=p.returncode,total_wall_seconds=time.monotonic()-start)
R.write_text(json.dumps(r,indent=2)+'\n')
print(out,end='');print(json.dumps(r));raise SystemExit(p.returncode)
