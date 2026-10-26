"""One owned child, CPU60, wall<=180 including cleanup and receipt."""
from pathlib import Path
import datetime, json, os, resource, signal, subprocess, sys, time
P=Path(__file__).resolve().parent
receipt=P/'run_receipt.json'
if receipt.exists() or (P/'results.json').exists(): raise SystemExit('one attempt only')
start=time.monotonic();cpu=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(TimeoutError('supervisor wall175')))
signal.alarm(175)
def limit():
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
child=None;status='UNSTARTED';out='';err='';code=None
try:
    child=subprocess.Popen([sys.executable,str(P/'check_fixed.py')],cwd=P,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True,preexec_fn=limit)
    try:
        out,err=child.communicate(timeout=160);code=child.returncode
        status='COMPLETED' if code==0 else 'FIRST_FAILURE_STOP'
    except (subprocess.TimeoutExpired,TimeoutError):
        status='TIMEOUT_NO_RETRY'
        os.killpg(child.pid,signal.SIGKILL)
        out,err=child.communicate(timeout=8);code=child.returncode
finally:
    if child is not None and child.poll() is None:
        os.killpg(child.pid,signal.SIGKILL);child.wait(timeout=5)
    signal.alarm(0)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    r={'status':status,'started_utc':started_utc,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'guard_pid':os.getpid(),'child_pid':None if child is None else child.pid,'exit_code':code,'timeout_seconds':160,'cpu_limit_seconds':60,'program_wall_hard_seconds':180,'cleanup_included':True,'guard_elapsed_seconds':time.monotonic()-start,'guard_cpu_seconds':time.process_time()-cpu,'child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'reaped':child is None or child.poll() is not None,'automatic_retries':0,'stdout':out,'stderr':err}
    receipt.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r))
if code!=0:raise SystemExit(1)
