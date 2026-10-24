#!/usr/bin/env python3
"""One predeclared discovery invocation; 300s hard wall cap, no retry."""
import hashlib,json,os,signal,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
out=ROOT/'discovery.run.json'
if out.exists() or (ROOT/'result.json').exists():
    raise SystemExit('Existing run/result; automatic repetition forbidden.')
source=ROOT/'search.py'
record={'command':[sys.executable,str(source)],'hard_seconds':300,
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'status':'started'}
out.write_text(json.dumps(record,indent=2)+'\n')
start=time.monotonic()
with (ROOT/'discovery.log').open('w') as stream:
    proc=subprocess.Popen(record['command'],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
    try:
        proc.wait(timeout=max(.001,300-(time.monotonic()-start)))
        record['status']='returned'
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid,signal.SIGKILL)
        proc.wait()
        record['status']='hard_timeout_killed'
# Reap only the process group this wrapper created, including a batch child.
try:
    os.killpg(proc.pid,0)
except ProcessLookupError:
    record['owned_group_cleanup']='already_absent'
else:
    os.killpg(proc.pid,signal.SIGKILL)
    record['owned_group_cleanup']='remaining_owned_group_killed'
record['source_sha256_after']=hashlib.sha256(source.read_bytes()).hexdigest()
record['source_unchanged']=record['source_sha256']==record['source_sha256_after']
record.update(exit_code=proc.returncode,wall_seconds=time.monotonic()-start,
              log_sha256=hashlib.sha256((ROOT/'discovery.log').read_bytes()).hexdigest())
out.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
