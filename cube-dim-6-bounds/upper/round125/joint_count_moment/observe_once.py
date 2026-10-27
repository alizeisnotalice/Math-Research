#!/usr/bin/env python3
"""Observe exactly one authorized125 child through full exit."""
import hashlib,json,os,resource,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
rel=json.loads((HERE/'parent_release.json').read_text())
assert rel.get('authorized')is True and str(rel['round'])=='125'
assert rel['observer_sha256']==sha(__file__) and rel['script_sha256']==sha(HERE/'run_prepared.py')
with(HERE/'observer_started.json').open('x')as f:
    json.dump({'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat()},f);f.write('\n')
before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
p=subprocess.run(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(HERE/'run_prepared.py')],cwd=str(HERE),capture_output=True,text=True)
elapsed=time.monotonic()-start;after=resource.getrusage(resource.RUSAGE_CHILDREN)
d={'exit_code':p.returncode,'complete_science_wall_seconds':elapsed,
   'complete_science_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
   'scope':'spawn before sole scientific child through final receipt, interpreter exit and reap; observer own write/exit excluded',
   'observer_pid':os.getpid(),'stdout':p.stdout,'stderr':p.stderr,'child_reaped':True,'children_created':1,
   'utc':datetime.now(timezone.utc).isoformat(),'wall_limit':None,'release_sha256':sha(HERE/'parent_release.json')}
t=HERE/'observed_tool_receipt.json.own-tmp'
with t.open('w')as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
os.replace(t,HERE/'observed_tool_receipt.json')
print(json.dumps({'exit_code':p.returncode,'complete_science_wall_seconds':elapsed,
                  'complete_science_CPU_seconds':d['complete_science_CPU_seconds'],
                  'observer_pid':os.getpid(),'receipt_sha256':sha(HERE/'observed_tool_receipt.json')}),flush=True)
raise SystemExit(p.returncode)
