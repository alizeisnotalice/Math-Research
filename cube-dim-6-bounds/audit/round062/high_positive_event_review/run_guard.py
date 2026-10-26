"""One independent audit child with an actual 180-second timeout, no restart."""
import datetime,json,os,signal,subprocess,sys,time
from pathlib import Path
P=Path(__file__).resolve().parent
assert not (P/'process_receipt.json').exists(),'one registered process only; no automatic retry'
command=[sys.executable,'-B',str(P/'review.py'),'--run']
r={'status':'STARTING','command':command,'wall_timeout_seconds':180,'automatic_retries':0,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
def save():(P/'process_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
save();tic=time.monotonic()
p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
r.update(status='RUNNING',owned_pid=p.pid,owned_process_group=p.pid);save()
try:
 out,err=p.communicate(timeout=180)
 r['status']='COMPLETED' if p.returncode==0 else 'FAILED_NO_RETRY'
except subprocess.TimeoutExpired:
 os.killpg(p.pid,signal.SIGTERM)
 try:out,err=p.communicate(timeout=5)
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGKILL);out,err=p.communicate();r['needed_SIGKILL']=True
 r['status']='TIMEOUT_STOP_NO_RETRY'
r.update(exit_code=p.returncode,stdout=out,stderr=err,elapsed_seconds=time.monotonic()-tic,child_reaped=True,active_handles=[],finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save();print(json.dumps(r,indent=2));raise SystemExit(0 if r['status']=='COMPLETED' else 2)
