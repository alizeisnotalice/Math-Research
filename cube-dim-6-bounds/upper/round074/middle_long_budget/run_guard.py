"""One owned child;180 seconds includes termination/reaping allowance."""
from pathlib import Path
import subprocess,time,json,os,sys
H=Path(__file__).resolve().parent
receipt=H/'process_receipt.json'
if receipt.exists():raise RuntimeError('Single invocation only')
start=time.monotonic()
p=subprocess.Popen([sys.executable,'-B',str(H/'evaluate.py')],cwd=H,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
r={'owned_pid':p.pid,'guard_pid':os.getpid(),'status':'RUNNING','wall_total_limit_seconds':180,'child_wait_limit_seconds':178,'automatic_retries':0}
receipt.write_text(json.dumps(r,indent=2)+'\n')
try:
    output,_=p.communicate(timeout=178)
except subprocess.TimeoutExpired:
    p.kill();output,_=p.communicate(timeout=2);r['timeout']=True
r.update(status='EXITED_AND_REAPED',exit_code=p.returncode,total_wall_seconds=time.monotonic()-start)
receipt.write_text(json.dumps(r,indent=2)+'\n')
print(output,end='');print(json.dumps(r));raise SystemExit(p.returncode)
