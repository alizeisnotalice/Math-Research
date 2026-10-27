#!/usr/bin/env python3
"""Single detached band-seed observer; logs and complete exit receipt."""
import hashlib,json,os,resource,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rel=json.loads((R/'parent_release.json').read_text())
assert rel.get('authorized')is True and str(rel.get('round'))=='137'
assert rel['observer_sha256']==sha(__file__)and rel['script_sha256']==sha(R/'run_prepared.py')
n=os.getpriority(os.PRIO_PROCESS,0)
if n<10:os.nice(10-n)
with(R/'observer_started.json').open('x')as f:json.dump({'pid':os.getpid(),'unix_time':time.time(),'release_sha256':sha(R/'parent_release.json')},f);f.write('\n');f.flush();os.fsync(f.fileno())
before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
with(R/'science_stdout.txt').open('x')as out,(R/'science_stderr.txt').open('x')as err:
 p=subprocess.run(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(R/'run_prepared.py')],cwd=R,stdout=out,stderr=err)
after=resource.getrusage(resource.RUSAGE_CHILDREN)
d={'exit_code':p.returncode,'complete_science_wall_seconds':time.monotonic()-start,'complete_science_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'scope':'spawn sole band-seed controller child through exit/reap; observer receipt excluded','child_reaped':True,'children_created':1,'priority':os.getpriority(os.PRIO_PROCESS,0),'time_or_CPU_ceiling':None,'automatic_restarts':0,'stdout_sha256':sha(R/'science_stdout.txt'),'stderr_sha256':sha(R/'science_stderr.txt'),'release_sha256':sha(R/'parent_release.json')}
t=R/'observed_tool_receipt.json.own-tmp'
with t.open('w')as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
os.replace(t,R/'observed_tool_receipt.json');print(json.dumps(d),flush=True);raise SystemExit(p.returncode)
