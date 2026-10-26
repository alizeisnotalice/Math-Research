"""One owned child, CPU60; wall175 with cleanup reserved inside180."""
from pathlib import Path
import datetime,json,os,resource,signal,subprocess,sys,time
P=Path(__file__).resolve().parent
receipt=P/'run_receipt.json'
if receipt.exists() or (P/'results.json').exists():raise SystemExit('one attempt only')
start=time.monotonic();cpu=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
started=datetime.datetime.now(datetime.timezone.utc).isoformat()
signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('guard wall175')))
signal.alarm(175)
def limits():resource.setrlimit(resource.RLIMIT_CPU,(60,60))
child=None;status='UNSTARTED';out='';err='';code=None
try:
 child=subprocess.Popen([sys.executable,str(P/'check_fixed.py')],cwd=P,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True,preexec_fn=limits)
 try:
  out,err=child.communicate(timeout=160);code=child.returncode;status='COMPLETED'if code==0 else'FIRST_FAILURE_STOP'
 except(subprocess.TimeoutExpired,TimeoutError):
  status='TIMEOUT_NO_RETRY';os.killpg(child.pid,signal.SIGKILL);out,err=child.communicate(timeout=8);code=child.returncode
finally:
 if child is not None and child.poll()is None:os.killpg(child.pid,signal.SIGKILL);child.wait(timeout=5)
 signal.alarm(0);after=resource.getrusage(resource.RUSAGE_CHILDREN)
 data={'status':status,'started_utc':started,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'guard_pid':os.getpid(),'child_pid':None if child is None else child.pid,'exit_code':code,'child_wait_seconds':160,'CPU_limit_seconds':60,'program_wall_limit_including_cleanup_receipt_seconds':180,'guard_seconds_after_child_cleanup_before_final_receipt_write':time.monotonic()-start,'guard_cpu_seconds':time.process_time()-cpu,'child_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'reaped':child is None or child.poll()is not None,'retries':0,'stdout':out,'stderr':err}
 receipt.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data))
if code!=0:raise SystemExit(1)
