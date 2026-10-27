#!/usr/bin/env python3
import sys,os,time,resource,runpy,traceback
from pathlib import Path
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import run_batch as R
if __name__=='__main__':
 t=time.monotonic();a=resource.getrusage(resource.RUSAGE_SELF);code=1;failure=None
 try:sys.argv=[str(P/'run_batch.py')];runpy.run_path(sys.argv[0],run_name='__main__');code=0
 except BaseException:failure=traceback.format_exc();print(failure,file=sys.stderr,flush=True)
 finally:
  b=resource.getrusage(resource.RUSAGE_SELF);R.atomic(P/'detached_exit_receipt.json',{'pid':os.getpid(),'exit_code':code,'failure':failure,'full_runner_wall_seconds':time.monotonic()-t,'self_CPU_user_seconds':b.ru_utime-a.ru_utime,'self_CPU_system_seconds':b.ru_stime-a.ru_stime,'native_calls':0,'exp_calls':0})
 raise SystemExit(code)
