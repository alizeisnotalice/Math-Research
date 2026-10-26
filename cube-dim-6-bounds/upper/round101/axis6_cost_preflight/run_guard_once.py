#!/usr/bin/env python3
"""Own-child guard. Default hard fuse stays armed through receipt/stdout/exit."""
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,177)
import os
import sys
import time
import resource
import subprocess
import json
from pathlib import Path

START=time.monotonic()
resource.setrlimit(resource.RLIMIT_CPU,(2,2))
signal.signal(signal.SIGXCPU,signal.SIG_DFL)
H=Path(__file__).resolve().parent
if (H/'run_receipt.json').exists() or (H/'results.json').exists():
    raise RuntimeError('exclusive once-only guard refuses existing result/receipt')
rec={'status':'RUNNING','guard_pid':os.getpid(),'outer_hard_fuse_seconds':177,
     'hard_fuse_disposition':'SIG_DFL; never cancelled',
     'CPU_guard_hard_seconds':2,'CPU_child_hard_seconds':28,'CPU_combined_hard_seconds':30,
     'child_hard_fuse_seconds':150,'child_wait_seconds':160,'own_group_cleanup_seconds':5,
     'total_wall_seconds_cap_including_final_io_exit':180,'automatic_retries':0}
p=subprocess.Popen([sys.executable,'-B',str(H/'check_exact.py')],cwd=H,
    stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
rec['owned_child_pid']=p.pid
rec['owned_child_pgid']=p.pid
stdout=stderr=''
try:
    stdout,stderr=p.communicate(timeout=160)
    rec.update(status='EXITED_AND_REAPED',child_reaped=True,exit_code=p.returncode)
except subprocess.TimeoutExpired:
    rec['timeout']=True
    try:os.killpg(p.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    try:
        stdout,stderr=p.communicate(timeout=5)
        rec.update(status='TIMED_OUT_OWN_GROUP_KILLED_REAPED',child_reaped=True,exit_code=p.returncode)
    except subprocess.TimeoutExpired:
        rec.update(status='FAILED_BOUNDED_CLEANUP',child_reaped=False,exit_code=1)
ru=resource.getrusage(resource.RUSAGE_CHILDREN)
rec.update(stdout=stdout,stderr=stderr,child_OS_CPU_seconds=ru.ru_utime+ru.ru_stime,
    guard_CPU_seconds=time.process_time(),wall_seconds_before_final_receipt_write=time.monotonic()-START,
    hard_fuse_still_armed=signal.getitimer(signal.ITIMER_REAL)[0]>0)
with (H/'run_receipt.json').open('x') as f:json.dump(rec,f,ensure_ascii=False,indent=2);f.write('\n')
print(json.dumps(rec,ensure_ascii=False),flush=True)
sys.stdout.flush();sys.stderr.flush()
# No timer cancellation, no unbounded atexit work. Tool wall is saved separately.
os._exit(0 if rec.get('child_reaped') and rec.get('exit_code')==0 else 1)
