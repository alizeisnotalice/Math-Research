#!/usr/bin/env python3
"""One observer child with direct detached logs; no runtime ceiling."""
import json,os,resource,subprocess,time
from pathlib import Path
from box_core import sha,atomic,require
HERE=Path(__file__).resolve().parent
rel=json.loads((HERE/'parent_release.json').read_text())
require(rel.get('authorized')is True and str(rel.get('round'))=='130','primary130 release')
require(rel['observer_sha256']==sha(__file__)and rel['script_sha256']==sha(HERE/'run_prepared.py'),'released observer/source')
priority=os.getpriority(os.PRIO_PROCESS,0)
if priority<10:os.nice(10-priority)
with(HERE/'observer_started.json').open('x')as f:json.dump({'pid':os.getpid(),'unix_time':time.time(),'release_sha256':sha(HERE/'parent_release.json')},f);f.write('\n');f.flush();os.fsync(f.fileno())
start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with(HERE/'science_stdout.txt').open('x')as out,(HERE/'science_stderr.txt').open('x')as err:
    result=subprocess.run(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(HERE/'run_prepared.py')],cwd=HERE,stdout=out,stderr=err)
after=resource.getrusage(resource.RUSAGE_CHILDREN)
receipt={'exit_code':result.returncode,'complete_science_wall_seconds':time.monotonic()-start,
         'complete_science_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
         'scope':'spawn sole science child through exit and reap; observer own receipt excluded','child_reaped':True,'children_created':1,
         'priority':os.getpriority(os.PRIO_PROCESS,0),'wall_or_CPU_limit':None,'automatic_restarts':0,
         'stdout_sha256':sha(HERE/'science_stdout.txt'),'stderr_sha256':sha(HERE/'science_stderr.txt'),'release_sha256':sha(HERE/'parent_release.json')}
atomic(HERE/'observed_tool_receipt.json',receipt);print(json.dumps(receipt),flush=True);raise SystemExit(result.returncode)
