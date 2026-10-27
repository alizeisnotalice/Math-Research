"""Terra-only once observer. Exact original-child wait4; no timeout/retry, complete lifecycle receipt."""
from pathlib import Path
import os,sys,time,json,subprocess,hashlib,traceback
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
import check_saved as C
def run():
 C.need(P==C.TARGET,'TERRA_ONLY_ROOT');release=C.read(P/'parent_release.json')
 for k,v in C.release_requirements().items():C.need(release.get(k)==v,'OBSERVER_RELEASE:'+k)
 C.need(release.get('main_approved')is True and release.get('independent_review_approved')is True,'OBSERVER_DUAL_GATE')
 C.need(not(P/'started.json').exists(),'NO_STARTED_AUDIT_OR_RESTART')
 C.save('observer_started.json',{'pid':os.getpid(),'release_sha256':C.sha(P/'parent_release.json'),'no_automatic_retry':True})
 os.nice(max(0,15-os.getpriority(os.PRIO_PROCESS,0)))
 for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[n]='1'
 start=time.monotonic();pid=None;code=None;usage=None;failure=None;reaped=False
 try:
  with(P/'checker_stdout.txt').open('xb')as out,(P/'checker_stderr.txt').open('xb')as err:
   child=subprocess.Popen([sys.executable,'-B',str(P/'check_saved.py')],stdin=subprocess.DEVNULL,stdout=out,stderr=err,cwd=P);pid=child.pid
   C.save('original_process.json',{'observer_pid':os.getpid(),'checker_pid':pid,'checker_sha256':C.sha(P/'check_saved.py'),'registration_sha256':C.sha(P/'registration.json')})
   waited,status,usage=os.wait4(pid,0);C.need(waited==pid,'EXACT_CHILD_REAP');code=os.waitstatus_to_exitcode(status);child.returncode=code;reaped=True
 except BaseException:
  failure=traceback.format_exc()
  # Preserve a still-owned original child; never dispatch a second one or kill it for time.
  if pid is not None and not reaped:
   waited,status,usage=os.wait4(pid,0);code=os.waitstatus_to_exitcode(status);reaped=True
 elapsed=time.monotonic()-start;ex=C.read(P/'exit_receipt.json')if(P/'exit_receipt.json').exists()else None
 valid=code==0 and reaped and failure is None and ex and ex.get('status')=='PASS_COMPLETE_SAVED156_INDEPENDENT_SCALAR768_AND_RATIONAL_CAPS'and(P/'certificate.json').exists()and C.sha(P/'certificate.json')==ex.get('certificate_sha256')and(P/'checker_stderr.txt').read_bytes()==b''
 receipt={'status':'TERMINAL_PASS'if valid else'FAILED_OR_INCOMPLETE_NO_RETRY','observer_pid':os.getpid(),'checker_pid':pid,'exit_code':code,'child_reaped':reaped,'full_child_wall_seconds':elapsed,'child_CPU_user_seconds':usage.ru_utime if usage else None,'child_CPU_system_seconds':usage.ru_stime if usage else None,'maxrss_platform_units':usage.ru_maxrss if usage else None,'maxrss_units':'bytes'if sys.platform=='darwin'else'KiB','failure':failure,'counts_final':ex['counts']if valid else None,'scalar_calls_final':ex['scalar_calls']if valid else None,'incomplete_count_policy':'uncheckpointed final work UNKNOWN; keep completed cone/scalar checkpoint snapshots','stdout_sha256':C.sha(P/'checker_stdout.txt')if(P/'checker_stdout.txt').exists()else None,'stderr_sha256':C.sha(P/'checker_stderr.txt')if(P/'checker_stderr.txt').exists()else None,'release_sha256':C.sha(P/'parent_release.json'),'automatic_retries':0,'producer_or_old_baseline_calls':0,'new_C6_bound':False}
 C.save('operator_terminal_receipt.json',receipt);print(json.dumps({'status':receipt['status'],'checker_pid':pid,'exit_code':code,'wall':elapsed}));return 0 if valid else 1
if __name__=='__main__':raise SystemExit(run())
