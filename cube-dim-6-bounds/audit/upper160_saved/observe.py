"""Terra-only single audit observer; no timeout/retry, exact controller wait4 including reaped workers."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time,traceback
sys.dont_write_bytecode=True
import check_saved as C
P=Path(__file__).resolve().parent
def main():
 C.need(P==C.TARGET,'TERRA_ROOT_ONLY');release=C.read(P/'parent_release.json')
 for k,v in C.release_requirements().items():C.need(release.get(k)==v,'OBSERVER_RELEASE:'+k)
 C.need(release.get('main_approved')is True and release.get('independent_review_approved')is True,'DUAL_RELEASE')
 C.need(not(P/'started.json').exists(),'NO_RESTART')
 C.save('observer_started.json',{'pid':os.getpid(),'release_sha256':C.sha(P/'parent_release.json'),'automatic_retries':0})
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 for k in C.ENV:os.environ[k]='1'
 start=time.monotonic();pid=None;rc=None;usage=None;failure=None;reaped=False;child=None
 try:
  with(P/'checker_stdout.txt').open('xb')as out,(P/'checker_stderr.txt').open('xb')as err:
   child=subprocess.Popen([sys.executable,'-B',str(P/'check_saved.py')],stdin=subprocess.DEVNULL,stdout=out,stderr=err,cwd=P);pid=child.pid
   C.save('original_process.json',{'observer_pid':os.getpid(),'checker_pid':pid,'checker_sha256':C.sha(P/'check_saved.py'),'registration_sha256':C.sha(P/'registration.json')})
   waited,status,usage=os.wait4(pid,0);C.need(waited==pid,'EXACT_ORIGINAL_CHECKER_WAIT4');rc=os.waitstatus_to_exitcode(status);child.returncode=rc;reaped=True
 except BaseException:
  failure=traceback.format_exc()
  if pid is not None and not reaped:
   try:
    waited,status,usage=os.wait4(pid,0);C.need(waited==pid,'EXACT_ORIGINAL_REAP_AFTER_ERROR');rc=os.waitstatus_to_exitcode(status);child.returncode=rc;reaped=True
   except BaseException:failure+='\nREAP_INCOMPLETE_UNKNOWN\n'+traceback.format_exc()
 elapsed=time.monotonic()-start;ex=C.read(P/'exit_receipt.json')if(P/'exit_receipt.json').exists()else None
 valid=rc==0 and reaped and failure is None and ex and ex.get('status')==C.PASS and ex.get('workers_issued')==ex.get('workers_reaped')==15 and (P/'certificate.json').exists()and C.sha(P/'certificate.json')==ex.get('certificate_sha256')and(P/'checker_stderr.txt').read_bytes()==b''
 receipt={'status':'TERMINAL_PASS'if valid else'FAILED_OR_INCOMPLETE_NO_RETRY','observer_pid':os.getpid(),'checker_pid':pid,'exit_code':rc,'child_reaped':reaped,'full_child_wall_seconds':elapsed,'wait4_user_CPU_seconds':usage.ru_utime if usage else None,'wait4_system_CPU_seconds':usage.ru_stime if usage else None,'CPU_scope':'controller plus workers actually waited by that original controller; no inference about unobserved children','failure':failure,'counts':ex['counts']if valid else None,'certificate_sha256':C.sha(P/'certificate.json')if valid else None,'checker_exit_sha256':C.sha(P/'exit_receipt.json')if ex else None,'stdout_sha256':C.sha(P/'checker_stdout.txt')if(P/'checker_stdout.txt').exists()else None,'stderr_sha256':C.sha(P/'checker_stderr.txt')if(P/'checker_stderr.txt').exists()else None,'release_sha256':C.sha(P/'parent_release.json'),'unobserved_work':'UNKNOWN if missing worker/exit/reap receipt; no retry','automatic_retries':0,'new_fifthroot_Arb_kernel_native_optimizer_old_baseline_calls':0,'new_C6_bound':False}
 C.save('observer_terminal_receipt.json',receipt);print(json.dumps({'status':receipt['status'],'checker_pid':pid,'exit_code':rc,'wall':elapsed}),flush=True)
 return 0 if valid else 1
if __name__=='__main__':raise SystemExit(main())
