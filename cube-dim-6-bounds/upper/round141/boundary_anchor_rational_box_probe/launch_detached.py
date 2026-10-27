#!/usr/bin/env python3
import os,sys,subprocess
from pathlib import Path
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import run_batch as R
if __name__=='__main__':
 R.release_gate();R.require(not(P/'started.json').exists(),'ALREADY_STARTED');os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));R.write_new(P/'detached_owner_lock.json',{'launcher_pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'automatic_restart':False})
 with(P/'detached.log').open('xb')as log:
  child=subprocess.Popen([sys.executable,str(P/'detached_runner.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=P,start_new_session=True,close_fds=True)
 R.write_new(P/'detached_launch_receipt.json',{'pid':child.pid,'start_new_session':True,'run_batch_sha256':R.sha(P/'run_batch.py'),'time_or_CPU_ceiling':None,'native_calls':0,'exp_calls':0});print(child.pid)
