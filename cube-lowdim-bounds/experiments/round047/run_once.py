"""One approved diagnostic only; external 300s hard process-group deadline."""
from pathlib import Path
from time import monotonic,time
from hashlib import sha256
import json,os,signal,subprocess,sys
root=Path(__file__).resolve().parent
meta=root/'search.run.json';log=root/'search.log';script=root/'search.py'
if any(p.exists() for p in (meta,log,root/'result.json',root/'solver_result.json')):
 raise SystemExit('Existing evidence: do not repeat or overwrite.')
assert sha256(script.read_bytes()).hexdigest()=='715678d73fa80e7eb108d9642a6b02fde242773be7feec8d44931359510731d2'
assert sha256((root/'solver.m').read_bytes()).hexdigest()=='d37a5f3194e8fce144f7eef94a2c2f84de03579180a68769abecf1c161ddcd8c'
start=monotonic();record={'status':'started','budget_seconds':300,'started_unix':time(),
 'command_relative':['python','-B','search.py'],'source_sha256_before':sha256(script.read_bytes()).hexdigest()}
meta.write_text(json.dumps(record,indent=2)+'\n')
with log.open('w') as stream:
 p=subprocess.Popen([sys.executable,'-B',str(script)],cwd=root,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
 record['pid']=p.pid;meta.write_text(json.dumps(record,indent=2)+'\n')
 try:
  code=p.wait(timeout=max(.001,300-(monotonic()-start)))
  record.update(status='terminal',exit_code=code)
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGKILL);p.wait()
  record.update(status='timed_out',exit_code=p.returncode)
try:os.killpg(p.pid,signal.SIGKILL)
except ProcessLookupError:pass
record.update(elapsed_seconds=monotonic()-start,log_sha256=sha256(log.read_bytes()).hexdigest(),
 source_sha256_after=sha256(script.read_bytes()).hexdigest())
meta.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
