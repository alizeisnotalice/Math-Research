from pathlib import Path
import subprocess,time,json,os,sys,signal
H=Path(__file__).resolve().parent;path=H/'process_receipt.json'
if path.exists():raise RuntimeError('Single invocation only')
start=time.monotonic();p=subprocess.Popen([sys.executable,'-B',str(H/'check_independent.py')],cwd=H,start_new_session=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
r={'owned_pid':p.pid,'guard_pid':os.getpid(),'status':'RUNNING','wall_total_limit_seconds':180,'work_timeout_seconds':170,'cleanup_reserve_seconds':10,'automatic_retries':0}
path.write_text(json.dumps(r,indent=2)+'\n')
try:output,_=p.communicate(timeout=170)
except subprocess.TimeoutExpired:
 os.killpg(p.pid,signal.SIGKILL);output,_=p.communicate(timeout=10);r['timeout']=True
r.update(status='EXITED_AND_REAPED',exit_code=p.returncode,total_wall_seconds=time.monotonic()-start)
path.write_text(json.dumps(r,indent=2)+'\n');print(output,end='');print(json.dumps(r));raise SystemExit(p.returncode)
