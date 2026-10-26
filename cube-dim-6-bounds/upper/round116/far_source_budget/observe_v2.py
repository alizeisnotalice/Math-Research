#!/usr/bin/env python3
"""Observe exactly one authorized scientific child through complete exit."""
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
PYTHON = '/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

release = json.loads((HERE/'parent_release_v2.json').read_text())
assert release['authorized'] is True and release['observer_sha256'] == sha(__file__)
assert release['script_sha256'] == sha(HERE/'resume_tail_v2.py')
with (HERE/'observer_started_v2.json').open('x') as f:
    json.dump({'pid': os.getpid(), 'utc': datetime.now(timezone.utc).isoformat(),
               'observer_sha256': sha(__file__)}, f, indent=2)
    f.write('\n')
before = resource.getrusage(resource.RUSAGE_CHILDREN)
start = time.monotonic()
p = subprocess.run([PYTHON, '-B', str(HERE/'resume_tail_v2.py')],
                   cwd=str(HERE), capture_output=True, text=True)
elapsed = time.monotonic()-start
after = resource.getrusage(resource.RUSAGE_CHILDREN)
data = {'scientific_child_exit_code': p.returncode,
        'child_complete_wall_seconds': elapsed,
        'child_CPU_seconds_complete': after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        'scope': 'before spawning sole child through child wait/reap after its final receipt and interpreter exit; observer own receipt write/exit excluded',
        'observer_pid': os.getpid(), 'utc': datetime.now(timezone.utc).isoformat(),
        'stdout': p.stdout, 'stderr': p.stderr, 'own_children_created': 1,
        'child_reaped': True, 'runtime_limit': None,
        'entry_sha256': sha(HERE/'resume_tail_v2.py'),
        'release_sha256': sha(HERE/'parent_release_v2.json')}
q = HERE/'observed_tool_receipt_v2.json.own-tmp'
with q.open('w') as f:
    json.dump(data, f, ensure_ascii=False, indent=2)
    f.write('\n')
    f.flush()
    os.fsync(f.fileno())
os.replace(q, HERE/'observed_tool_receipt_v2.json')
print(json.dumps(data, ensure_ascii=False), flush=True)
raise SystemExit(p.returncode)
