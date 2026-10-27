#!/usr/bin/env python3
import os,sys,time,resource,runpy,traceback,datetime
from pathlib import Path
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import run_batch as R
if __name__=='__main__':
 start=time.monotonic();cpu=resource.getrusage(resource.RUSAGE_SELF);code=1;failure=None
 try:
  sys.argv=[str(P/'run_batch.py')];runpy.run_path(sys.argv[0],run_name='__main__');code=0
 except BaseException:failure=traceback.format_exc();print(failure,file=sys.stderr,flush=True)
 finally:
  end=resource.getrusage(resource.RUSAGE_SELF);R.atomic(P/'detached_exit_receipt.json',{'pid':os.getpid(),'exit_code':code,'failure':failure,'full_runner_wall_seconds':time.monotonic()-start,'self_cpu_user_seconds':end.ru_utime-cpu.ru_utime,'self_cpu_system_seconds':end.ru_stime-cpu.ru_stime,'cpu_scope':'same detached producer process; RUSAGE_SELF start/end snapshots, includes preflight and scientific body','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 raise SystemExit(code)
