#!/usr/bin/env python3
"""One authorized scientific child, observed through receipt and full exit."""
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time
from datetime import datetime, timezone
HERE=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
release=json.loads((HERE/'parent_release.json').read_text())
assert release['authorized'] is True and str(release['round'])=='120'
assert release['observer_sha256']==sha(__file__)
assert release['script_sha256']==sha(HERE/'run_prepared.py')
with (HERE/'observer_started.json').open('x') as f:
    json.dump({'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat()},f)
    f.write('\n')
ru0=resource.getrusage(resource.RUSAGE_CHILDREN)
start=time.monotonic()
p=subprocess.run(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(HERE/'run_prepared.py')],cwd=str(HERE),capture_output=True,text=True)
wall=time.monotonic()-start
ru1=resource.getrusage(resource.RUSAGE_CHILDREN)
data={'exit_code':p.returncode,'complete_science_wall_seconds':wall,
      'complete_science_CPU_seconds':ru1.ru_utime+ru1.ru_stime-ru0.ru_utime-ru0.ru_stime,
      'scope':'spawn before scientific child through final child receipt, interpreter exit and reap; observer own write/exit excluded',
      'observer_pid':os.getpid(),'stdout':p.stdout,'stderr':p.stderr,'child_reaped':True,
      'utc':datetime.now(timezone.utc).isoformat(),'children_created':1,'runtime_limit':None,
      'science_entry_sha256':sha(HERE/'run_prepared.py'),'release_sha256':sha(HERE/'parent_release.json')}
tmp=HERE/'observed_tool_receipt.json.own-tmp'
with tmp.open('w') as f:
    json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
os.replace(tmp,HERE/'observed_tool_receipt.json')
print(json.dumps(data,ensure_ascii=False),flush=True)
raise SystemExit(p.returncode)
