#!/usr/bin/env python3
import hashlib,json,os,resource,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
reg=json.loads((R/'registration.json').read_text());assert reg['observer_sha256']==sha(__file__)and reg['checker_sha256']==sha(R/'check_saved_independent.py')
with(R/'observer_started.json').open('x')as f:json.dump({'pid':os.getpid(),'unix_time':time.time(),'registration_sha256':sha(R/'registration.json')},f);f.write('\n');f.flush();os.fsync(f.fileno())
n=os.getpriority(os.PRIO_PROCESS,0)
if n<10:os.nice(10-n)
before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
p=subprocess.run(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(R/'check_saved_independent.py')],cwd=R,capture_output=True,text=True)
after=resource.getrusage(resource.RUSAGE_CHILDREN)
d={'exit_code':p.returncode,'complete_wall_seconds':time.monotonic()-start,'complete_child_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'stdout':p.stdout,'stderr':p.stderr,'child_reaped':True,'children_created':1,'scope':'spawn sole independent saved141 checker through exit/reap; observer own write excluded','priority':os.getpriority(os.PRIO_PROCESS,0),'checker_sha256':sha(R/'check_saved_independent.py'),'registration_sha256':sha(R/'registration.json'),'automatic_retries':0,'time_or_CPU_ceiling':None}
t=R/'actual_tool_receipt.json.own-tmp'
with t.open('w')as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
os.replace(t,R/'actual_tool_receipt.json');print(json.dumps(d),flush=True);raise SystemExit(p.returncode)
